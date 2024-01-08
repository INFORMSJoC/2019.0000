import torch
from RIM import RIMCell
import torch.nn as nn
from numpy import random
import torch.nn.functional as F
import torch.nn.init as init

# CUDA for PyTorch
use_cuda = torch.cuda.is_available()
device = torch.device("cuda:0" if use_cuda else "cpu")


def init_weight(n):
    if isinstance(n, nn.Linear):
        nn.init.xavier_uniform_(n.weight)
        nn.init.constant_(n.bias, 0)


def init_lstm(m):
    if isinstance(m, nn.LSTM):
        for param in m.parameters():
            if len(param.shape) >= 2:
                init.orthogonal_(param.data)
            else:
                init.normal_(param.data)


def linear(in_dim, out_dim, bias=True):
    lin = nn.Linear(in_dim, out_dim, bias=bias)
    nn.init.xavier_uniform_(lin.weight)
    if bias:
        lin.bias.data.zero_()

    return lin


class PerturbedTopK(nn.Module):
    def __init__(self, k: int, num_samples: int = 250, sigma: float = 0.05):
         super(PerturbedTopK, self).__init__()
         self.num_samples = num_samples
         self.sigma = sigma
         self.k = k

    def __call__(self, x):
        return PerturbedTopKFunction.apply(x, self.k, self.num_samples, self.sigma)


class PerturbedTopKFunction(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, k: int, num_samples: int = 1000, sigma: float = 0.05):
        b, d = x.shape
        # for Gaussian: noise and gradient are the same.
        noise = torch.normal(mean=0.0, std=1.0, size=(b, num_samples, d)).to(x.device)

        perturbed_x = x[:, None, :] + noise * sigma # b, nS, d
        topk_results = torch.topk(perturbed_x, k=k, dim=-1, sorted=False)
        indices = topk_results.indices # b, nS, k
        indices = torch.sort(indices, dim=-1).values # b, nS, k

        # b, nS, k, d
        perturbed_output = torch.nn.functional.one_hot(indices, num_classes=d).float()
        indicators = perturbed_output.mean(dim=1) # b, k, d

        # constants for backward
        ctx.k = k
        ctx.num_samples = num_samples
        ctx.sigma = sigma

        # tensors for backward
        ctx.perturbed_output = perturbed_output
        ctx.noise = noise

        return indicators

    @staticmethod
    def backward(ctx, grad_output):
        if grad_output is None:
         return tuple([None] * 5)

        noise_gradient = ctx.noise
        expected_gradient = (
                torch.einsum("bnkd,bnd->bkd", ctx.perturbed_output, noise_gradient)
                / ctx.num_samples
                / ctx.sigma
         )
        grad_input = torch.einsum("bkd,bkd->bd", grad_output, expected_gradient)
        return (grad_input,) + tuple([None] * 5)


class Retrieval_RIM(nn.Module):

    def __init__(self, input_dim, batch_size, hidden_dim, output_dim, num_layers, fc_dim, mode,
                 dropout_rate, topk, query_num):
        super(Retrieval_RIM, self).__init__()
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.fc_dim = fc_dim
        self.num_layers = num_layers
        self.mode = mode
        self.batch_size = batch_size
        self.lstm4heads = []
        self.query_num = query_num
        self.topk = topk

        self.token_self_attention = nn.Sequential(
            nn.Linear(self.input_dim, self.hidden_dim*self.num_layers),
            nn.Tanh(),
            nn.Linear(self.hidden_dim*self.num_layers, 1),
            nn.Tanh(),
        )

        self.self_attention = nn.Sequential(
            nn.Linear(self.input_dim, self.hidden_dim*self.num_layers),
            nn.Tanh(),
            nn.Linear(self.hidden_dim*self.num_layers, 1),
            nn.Tanh(),
        )
        # initialize attention
        init_weight(self.token_self_attention)
        init_weight(self.self_attention)

        RIMs = RIMCell(device, input_dim, hidden_size=self.hidden_dim, num_units=self.num_layers)

        self.bilinear_matrix = nn.ParameterList()

        self.lstm4heads = []

        # Three heads share the same RIMs parameters
        for heads in range(self.query_num):
            bilinear = torch.nn.Parameter(torch.randn(self.input_dim, self.input_dim))
            self.bilinear_matrix.append(bilinear)
            self.lstm4heads.append(RIMs)

        self.lstm4heads = torch.nn.ModuleList(self.lstm4heads)

        # RIMs for short-term memory
        self.lstm4today = RIMCell(device, input_dim, hidden_size=self.hidden_dim, num_units=self.num_layers)
        self.dropout = nn.Dropout(dropout_rate)

        self.fc = nn.Linear(self.hidden_dim * self.num_layers * (self.query_num + 1), self.fc_dim)
        self.tanh = nn.Tanh()

        # Define the output layer
        self.linear = nn.Linear(self.fc_dim, output_dim)

    def return_idx(self, x, k):
        topk_x = torch.topk(x, k=k, dim=-1, sorted=False)
        x_indices = topk_x.indices
        x_indices = torch.sort(x_indices, dim=-1).values
        return x_indices

    def init_rim(self, batch_size):
        return torch.randn(batch_size, self.num_layers, self.hidden_dim).cuda(), \
               torch.randn(batch_size, self.num_layers, self.hidden_dim).cuda()

    def Retrive_Mem(self, mem, one_day, batch_size, num_samples, sigma, ltm_event_mask, stm_event_mask):

        #  B->Batch, W-> sequence length, N-> number of news in window, H-> the dimension of input
        # (B * N * H) -> (B * N * 1) = [B, 20, 768] -> [B, 20, 1]
        weights = self.self_attention(one_day).squeeze()

        if stm_event_mask is not None:
            weights = weights.masked_fill(stm_event_mask, -9e15)

        query_w = F.softmax(weights, dim=1)

        # output perturb_topk
        short_term_topk = PerturbedTopKFunction.apply(query_w, self.query_num, num_samples, sigma)#self.stm_topk(query_w)
        short_term_topk_indices = self.return_idx(query_w, self.query_num)
        new_query = torch.matmul(short_term_topk, one_day)
        retrieved_vec = torch.tensor((), dtype=torch.float, requires_grad=False).to(device)
        lm_topk_indices = torch.tensor((), dtype=torch.float, requires_grad=False).to(device)

        for i in range(self.query_num):

            biMatrix = self.bilinear_matrix[i]
            queryMatrix = new_query[:, i, :].unsqueeze(-1)
            sim_mtx = torch.matmul(mem, torch.matmul(biMatrix.cuda(), queryMatrix.cuda())).squeeze()
            
            # mask empty events in ltm
            if ltm_event_mask is not None:
                sim_mtx = sim_mtx.masked_fill(ltm_event_mask, -9e15)

            w = F.softmax(sim_mtx, dim=1)
            lt_topk = PerturbedTopKFunction.apply(w, self.topk, num_samples, sigma)
            indices = self.return_idx(w, self.topk)
            lm_topk_indices = torch.cat((lm_topk_indices, indices.unsqueeze(1)), dim=1)
            new_mem = torch.matmul(lt_topk, mem)
            new_mem = torch.cat((new_mem, queryMatrix.squeeze().unsqueeze(1)), dim=1)

            hs, cs = self.init_rim(batch_size)
            xs = torch.split(new_mem, 1, 1)
            for x in xs:
                hs, cs = self.lstm4heads[i](x, hs, cs)

            lstm_output_onehead = hs.contiguous().view(batch_size, -1)
            retrieved_vec = torch.cat((retrieved_vec, lstm_output_onehead.to(device)), dim=1)

        return retrieved_vec, short_term_topk_indices, lm_topk_indices

    def forward(self, input, mem, num_samples, sigma, ltm_event_mask, stm_event_mask, ltm_token_mask, stm_token_mask):

        # [batch, days=30, titleLen=20, dim=768]->[batch, days, dim]
        input = input.squeeze()
        batch_size = input.size(0)
        # word_att size [batch, days, title_len, 1]
        word_att = self.token_self_attention(input)
        word_att = word_att.squeeze()

        word_wight = F.softmax(word_att, dim=2)
        short_sent_embd = (input * word_wight.unsqueeze(-1)).sum(dim=2)

        # # [batch, days=900, titleLen=20, dim=768]->[batch, days, dim]
        word_att_mem = self.token_self_attention(mem)
        word_att_mem = word_att_mem.squeeze()

        word_wight_mem = F.softmax(word_att_mem, dim=2)
        mem_sent_embd = (mem * word_wight_mem.unsqueeze(-1)).sum(dim=2)

        retrived_vec, st_topk_indices, lt_topk_indices = self.Retrive_Mem(mem_sent_embd, short_sent_embd, batch_size, num_samples, sigma, ltm_event_mask, stm_event_mask)

        hs, cs = self.init_rim(batch_size)
        xs = torch.split(short_sent_embd, 1, 1)
        for x in xs:
            hs, cs = self.lstm4today(x, hs, cs)

        today_lstm_output = hs.contiguous().view(batch_size, -1)
        lstm_outputs_collections = torch.cat((retrived_vec, today_lstm_output), dim=1)
        y_pred = self.linear(self.tanh(self.fc(self.dropout(lstm_outputs_collections))))

        return y_pred.view(-1), word_wight, word_wight_mem, st_topk_indices, lt_topk_indices


if __name__ == '__main__':
    print(torch.device('cuda'))
    model = Retrieval_RIM(768, batch_size=32, hidden_dim=64, output_dim=1,
                          num_layers=4, fc_dim=200, mode='train', dropout_rate=0.2).cuda()
    x = torch.randn(16, 30, 20, 768).cuda()
    mem = torch.randn(16, 100, 20, 768).cuda()
    y = torch.randn(16, 1).cuda()
    model(x, mem)