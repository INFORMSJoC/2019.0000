import json
import os
import pickle
import shutil
import numpy as np
import pandas
import torch
from torch.utils.data.dataset import Dataset
from torch.utils import data

# CUDA for PyTorch
use_cuda = torch.cuda.is_available()
device = torch.device("cuda:0" if use_cuda else "cpu")


class Dataset(Dataset):
    # Characterizes a dataset for PyTorch
    def __init__(self, list_IDs):
        """Initialization"""
        self.list_IDs = list_IDs

    def __len__(self):
        # Denotes the total number of samples
        return len(self.list_IDs)

    def __getitem__(self, index):
        # Generates one sample of data
        memory = torch.tensor(()).to(device)
        memory_token = []
        SPR1_batch = torch.tensor(())
        SPR22_batch = torch.tensor(())
        SPVt1_batch = torch.tensor(())
        SPVt22_batch = torch.tensor(())
        shortMem = torch.tensor(()).to(device)
        shortMem_token = []
        ltm_text = []
        stm_text = []

        for ID in index:
            shortM, longM, SPR1, SPR22, SPVt1, SPVt22, shortTokens, longTokens, raw_ltm_text, raw_stm_text = torch.load(ID)

            # Long memory size is [batch, 900, 20, 768]
            memory = torch.cat((memory, longM.unsqueeze(0)), 0)

            # The long memory token size is [batch, 900, 20]
            memory_token.append(longTokens)
            SPR1_batch = torch.cat((SPR1_batch, SPR1.unsqueeze(0)), 0)
            SPR22_batch = torch.cat((SPR22_batch, SPR22.unsqueeze(0)), 0)
            SPVt1_batch = torch.cat((SPVt1_batch, SPVt1.unsqueeze(0)), 0)
            SPVt22_batch = torch.cat((SPVt22_batch, SPVt22.unsqueeze(0)), 0)

            # The short term memory
            shortMem = torch.cat((shortMem, shortM.unsqueeze(0)), dim=0)
            shortMem_token.append(shortTokens)

            ltm_text.append(raw_ltm_text)
            stm_text.append(raw_stm_text)

            del shortM
            del longM

        return shortMem, memory, SPR1_batch, SPR22_batch, SPVt1_batch, SPVt22_batch, shortMem_token, memory_token, ltm_text, stm_text

