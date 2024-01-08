import os
import json
import random
import pickle
import shutil
import argparse
import subprocess
import numpy as np
from model import *
from tqdm import tqdm
from torch.autograd import Variable as var
from WSJ_DataLoader import Dataset
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_fscore_support


def seed_torch(seed):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed) # if you are using multi-GPU.
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True


def event_mask(ltm_text, stm_text, ltm_token, stm_token):
    '''
    :param ltm_text: events in longterm memory
    :param stm_text: events in shorterm memory
    :param ltm_token: tokens in longterm memory
    :param stm_token: tokens in shorterm memory
    :return: ltm event mask, stm event mask, ltm token mask, stm token mask
    '''

    ltm_event_pad = np.array(ltm_text) == 'pad'
    stm_event_pad = np.array(stm_text) == 'pad'
    ltm_event_pad = torch.tensor(ltm_event_pad).to(device)
    stm_event_pad = torch.tensor(stm_event_pad).to(device)

    ltm_token_pad = np.array(ltm_token) == 'pad'
    stm_token_pad = np.array(stm_token) == 'pad'
    ltm_token_pad = torch.tensor(ltm_token_pad).to(device)
    stm_token_pad = torch.tensor(stm_token_pad).to(device)

    return ltm_event_pad, stm_event_pad, ltm_token_pad, stm_token_pad


def MAE(y_pred_train, y_true):
    return np.mean(np.abs(y_pred_train - y_true))


def MSE(y_pred_train, y_true):
    return np.mean((y_pred_train - y_true)**2)


def train(num_epochs, hidden_dim, batch_size, output_dim, num_layers, fc_dim, learning_rate, mode, lowercased,
          model_name, DAYS_BACK, time_window, dropout_rate, top_n, task, round, filter, data_size, seed, headnum):

    seed_torch(seed)
    # The output size of BERT
    lstm_input_size = 768

    if model_name == 'LASR-RIM':
        model = Retrieval_RIM(lstm_input_size, batch_size=batch_size, hidden_dim=hidden_dim, output_dim=output_dim,
                        num_layers=num_layers, fc_dim=fc_dim, mode=mode, dropout_rate=dropout_rate,  topk=top_n, query_num=headnum)

    loss_fn = torch.nn.MSELoss()
    optimiser = torch.optim.Adam(model.parameters(), lr=learning_rate)

    # CUDA for PyTorch
    use_cuda = torch.cuda.is_available()
    device = torch.device("cuda:0" if use_cuda else "cpu")

    if use_cuda:
        model = model.to(device)

    # load the index of datset
    train_data_path = '../train_dataset_30_LASER'
    print(train_data_path)
    train_data_list = sorted(map(np.long, os.listdir(train_data_path)))
    train_data_list_with_path = [train_data_path + '/'+ str(id) for id in train_data_list]

    val_data_path = '../val_dataset_30_LASER'
    val_data_list = sorted(map(np.long, os.listdir(val_data_path)))
    val_data_list_with_path = [val_data_path + '/' + str(id) for id in val_data_list]

    test_data_path = '../test_dataset_30_LASER'
    test_data_list = sorted(map(np.long, os.listdir(test_data_path)))
    test_data_list_with_path = [test_data_path + '/' + str(id) for id in test_data_list]

    # Generators
    training_set = Dataset(train_data_list_with_path)
    #training_generator = [train_data_list[i: i + batch_size] for i in range(0, len(train_data_list), batch_size)]

    validation_set = Dataset(val_data_list_with_path)
    validation_generator = [val_data_list_with_path[i: i + batch_size] for i in range(0, len(val_data_list_with_path), batch_size)]

    test_set = Dataset(test_data_list_with_path)
    test_generator = [test_data_list_with_path[i: i + batch_size] for i in range(0, len(test_data_list_with_path), batch_size)]

    num_samples = 250
    sigma = sigma_r = 0.5

    best_MAE_val = np.inf  # minimum value on val dataset
    patience = 30
    patience_count = 0

    mae_test = []
    mse_test = []

    pred_list_val = []
    true_list_val = []

    pred_list_test = []
    true_list_test = []

    for epoch in tqdm(range(num_epochs)):
        run_test = False
        np.random.shuffle(train_data_list_with_path)
        training_generator = [train_data_list_with_path[i: i + batch_size] for i in range(0, len(train_data_list_with_path), batch_size)]

        model.train()

        print(mode + '  Epoch ' + str(epoch) + ' start......')

        # Clear stored gradient
        model.zero_grad()

        for batch in training_generator:
            # Transfer to GPU
            local_batch_train, memory_train, SPR1_train, SPR22_train, SPVt1_train, SPVt22_train, shortMem_token_train, memory_token_train, ltm_text_train, stm_text_train = training_set.__getitem__(batch)

            if task == 'SPVt1_future':
                local_labels_train = SPVt1_train
            elif task == 'SPVt22_future':
                local_labels_train = SPVt22_train
            else:
                print('No legal dataset name!')

            local_labels_train = local_labels_train.squeeze()
            batch_s_train, batch_mem_train, batch_value_train = local_batch_train.to(device), memory_train.to(device), local_labels_train.to(device)

            ltm_event_mask_train, stm_event_mask_train, ltm_token_mask_train, stm_token_mask_train \
                 = event_mask(ltm_text_train, stm_text_train, memory_token_train, shortMem_token_train)

            y_pred_train, _, _, _, _ = model(batch_s_train, batch_mem_train, num_samples, sigma, ltm_event_mask_train, stm_event_mask_train, ltm_token_mask_train, stm_token_mask_train)

            loss = loss_fn(y_pred_train, batch_value_train)

            loss.backward()
            optimiser.step()
            optimiser.zero_grad()

        # Validation
        with torch.set_grad_enabled(False):
            print('starting validation')
            model.eval()

            for batch_val in validation_generator:
                # Transfer to GPU
                local_batch_val, memory_val, SPR1_val, SPR22_val, SPVt1_val, SPVt22_val, shortMem_token_val, memory_token_val, ltm_text_val, stm_text_val = validation_set.__getitem__(batch_val)

                if task == 'SPVt1_future':
                    local_labels_val = SPVt1_val
                elif task == 'SPVt22_future':
                    local_labels_val = SPVt22_val
                else:
                    print('No legal dataset name!')

                local_labels_val = local_labels_val.squeeze()
                batch_s_val, batch_memory_val, batch_value_val = local_batch_val.to(device), memory_val.to(device), local_labels_val.to(device)

                ltm_event_mask_val, stm_event_mask_val, ltm_token_mask_val, stm_token_mask_val \
                    = event_mask(ltm_text_val, stm_text_val, memory_token_val, shortMem_token_val)

                y_pred_val, _, _, _, _ = model(batch_s_val, batch_memory_val, 1, 0, ltm_event_mask_val, stm_event_mask_val, ltm_token_mask_val, stm_token_mask_val)

                pred_list_val += y_pred_val.float().tolist()
                true_list_val += batch_value_val.tolist()

            # map_list for test
            mae_val_value = MAE(np.array(pred_list_val), np.array(true_list_val))

            if mae_val_value < best_MAE_val:
                print('new best epoch.......')
                best_MAE_val = mae_val_value
                run_test = True
                patience_count = 0
            else:
                patience_count += 1

        if patience_count > patience:
            print("Early stopping criterion satisfied.")
            break

        # Testing
        if run_test:
            with torch.set_grad_enabled(False):
                print('starting testing')
                model.eval()
                test_explain_data = []
                for batch_test in test_generator:
                    # Transfer to GPU
                    local_batch_test, memory_test, SPR1_test, SPR22_test, SPVt1_test, SPVt22_test, shortMem_token_test, memory_token_test, ltm_text_test, stm_text_test = test_set.__getitem__(batch_test)

                    if task == 'SPVt1_future':
                        local_labels_test = SPVt1_test
                    elif task == 'SPVt22_future':
                        local_labels_test = SPVt22_test
                    else:
                        print('No legal dataset name!')

                    local_labels_test = local_labels_test.squeeze()
                    batch_s_test, batch_memory_test, batch_value_test = local_batch_test.to(device), memory_test.to(device), local_labels_test.to(device)

                    ltm_event_mask_test, stm_event_mask_test, ltm_token_mask_test, stm_token_mask_test \
                        = event_mask(ltm_text_test, stm_text_test, memory_token_test, shortMem_token_test)

                    y_pred_test, _, _, _, _ = model(batch_s_test, batch_memory_test, 1, 0, ltm_event_mask_test, stm_event_mask_test, ltm_token_mask_test, stm_token_mask_test)
                    pred_list_test += y_pred_test.float().tolist()
                    true_list_test += batch_value_test.tolist()

                mae_test.append(MAE(np.array(pred_list_test), np.array(true_list_test)))
                mse_test.append(MSE(np.array(pred_list_test), np.array(true_list_test)))

    print('Training and validation are completed!')


if __name__ == "__main__":
    # Create the parser and add arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='task', help="The task could be SPVt1_future or SPVt22_future")

    # Parse and print the results
    args = parser.parse_args()
    print(args.task)

    train(num_epochs=50, hidden_dim=75, batch_size=8, output_dim=1, num_layers=4, fc_dim=300,
          learning_rate=1e-4, mode='train', lowercased=True, model_name='LASR-RIM', DAYS_BACK=1,
          time_window=24, dropout_rate=0.2, top_n=30, task=args.task, round=0, filter=filter,
          data_size='large', seed=25535, headnum=3)