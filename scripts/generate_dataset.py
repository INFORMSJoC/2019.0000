import json
import os
import pickle
import shutil
import numpy as np
import pandas
import torch
from torch.utils.data.dataset import Dataset
from torch.utils import data
from pytorch_pretrained_bert import BertTokenizer, BertModel
# OPTIONAL: if you want to have more information on what's happening, activate the logger as follows
import logging
from tqdm import tqdm

# CUDA for PyTorch
use_cuda = torch.cuda.is_available()
# #use_cuda = False
device = torch.device("cuda:0" if use_cuda else "cpu")

logging.basicConfig(level=logging.INFO)

# Load pre-trained model tokenizer (vocabulary)
tokenizer = BertTokenizer.from_pretrained('bert-base-uncased', do_lower_case=True)
model = BertModel.from_pretrained('bert-base-uncased')
model.to(device)
model.eval()


class Dataset2(Dataset):
    # Characterizes a dataset for PyTorch
    def __init__(self, list_IDs, dataset_path):
        """Initialization"""
        self.list_IDs = list_IDs
        self.dataset_path = dataset_path

    def __len__(self):
        # Denotes the total number of samples
        return len(self.list_IDs)

    def convert_bert_embed(self, title):
        if title == 'pad':
            return torch.zeros(20, 768).cuda(), ['pad']*20
        else:
            processed_title = tokenizer.tokenize(title)
            tokenized_text = ["[CLS]"] + processed_title + ["[SEP]"]
            # Define sentence A and B indices associated to 1st and 2nd sentences (see paper)
            segments_ids = [0] * len(tokenized_text)
            indexed_tokens = tokenizer.convert_tokens_to_ids(tokenized_text)
            # Convert inputs to PyTorch tensors
            tokens_tensor = torch.tensor([indexed_tokens]).to(device)
            segments_tensors = torch.tensor([segments_ids]).to(device)
            input_mask = torch.tensor([1] * len(tokens_tensor)).to(device)
            # Get the pooled_output which could be considered as encoded sentence
            with torch.no_grad():
                all_layer, _ = model(tokens_tensor, segments_tensors, input_mask)
            last_layer = all_layer[-2].squeeze()[1:-1]
            title_token = tokenized_text[1:-1]
            sent_len = len(last_layer)
            if sent_len > 20:
                last_layer = last_layer[0:20]
                title_token = title_token[0:20]
            else:
                for i in range(20 - sent_len):
                    last_layer = torch.cat((last_layer, torch.zeros(1, 768).cuda()), dim=0)
                    title_token.append('pad')

        return last_layer, title_token

    def __getitem__(self, index):
        # Generates one sample of data

        longM, SPR1, SPR22, SPVt1, SPVt22, shortM = torch.load(self.dataset_path + '/' + str(index))
        each_memory = torch.tensor(()).to(device)
        each_token = []
        # convert M to size [900, 20, 768]
        for title in longM:
            # the size of each_title is [20, 768]
            each_title, tokens = self.convert_bert_embed(title)
            each_memory = torch.cat((each_memory, each_title.unsqueeze(0)), 0)
            each_token.append(tokens)

        each_shortM = torch.tensor(()).to(device)
        each_Stoken = []
        for Stitle in shortM:
            each_Stitle, Stoken = self.convert_bert_embed(Stitle)
            each_shortM = torch.cat((each_shortM, each_Stitle.unsqueeze(0)), dim=0)
            each_Stoken.append(Stoken)

        return (each_shortM, each_memory, SPR1, SPR22, SPVt1, SPVt22, each_Stoken, each_token, longM, shortM)
 

def clean_dir(data_path):
    if not os.path.exists(data_path):
        os.makedirs(data_path)
    else:
        shutil.rmtree(data_path)
        os.mkdir(data_path)


def initialize_dataset(time_window, top_n, data_size):

    # load the index of datset, E:/Pycharm workspace
    #mnt / Jie_workspace / LSMR / Multi - task - Generate - Dataset / Multi - task - Dataset - longterm - titles
    train_data_path = './train_dataset_' + str(time_window) + '_' + str(top_n) + '_' + str(data_size) + '_' + 'back_30'
    train_data_list = sorted(map(np.long, os.listdir(train_data_path)))
    training_set = Dataset2(train_data_list, train_data_path)

    val_data_path = './val_dataset_' + str(time_window) + '_' + str(top_n) + '_' + str(data_size) + '_' + 'back_30'
    val_data_list = sorted(map(np.long, os.listdir(val_data_path)))
    validation_set = Dataset2(val_data_list, val_data_path)
    
    test_data_path = './test_dataset_' + str(time_window) + '_' + str(top_n) + '_' + str(data_size) + '_' + 'back_30'
    test_data_list = sorted(map(np.long, os.listdir(test_data_path)))
    test_set = Dataset2(test_data_list, test_data_path)

    train_data_path = '../train_dataset_30_LASER'

    # clean previous data
    clean_dir(train_data_path)

    print('generating trainging dataset............................')
    for i in tqdm(range(30, 2898, 1)): #range(30, 2898, 1) #range(training_set.__len__())
        # if i < 10:
        #     continue
        torch.save(training_set.__getitem__(i), train_data_path + '/' + str(i))
    print('complete generating trainging dataset............................')

    print('generating val dataset............................')
    
    val_data_path = '../val_dataset_30_LASER'
    clean_dir(val_data_path)
    for i in tqdm(range(validation_set.__len__())):
        torch.save(validation_set.__getitem__(i), val_data_path + '/' + str(i))
    print('complete generating val dataset............................')
    
    print('generating test dataset............................')
    
    test_data_path = '../test_dataset_30_LASER'
    clean_dir(test_data_path)
    for i in tqdm(range(test_set.__len__())):
        torch.save(test_set.__getitem__(i), test_data_path + '/' + str(i))
    print('complete generating test dataset............................')


if __name__ == '__main__':

    initialize_dataset(time_window=24, top_n=30,  data_size='large')