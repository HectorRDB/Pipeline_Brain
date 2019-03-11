#!/opt/intel/intelpython3/bin/python3.6
import urllib.request
import os
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from scvi.models import *
from scvi.inference import UnsupervisedTrainer
from scvi.dataset import LoomDataset, CsvDataset

# Set the parameters
loc = "/pylon5/ib5phhp/hectorrb/ProcessedData/"
out = "/pylon5/ib5phhp/hectorrb/ProcessedData/
show_plot = True
n_epochs=400
lr=1e-3
use_batches=True

# Load the data
data = CsvDataset(filename = "10x_cells_MOp_filt.csv.gz",
                  save_path = loc,
                  compression = 'gzip',
                  batch_ids_file = "10x_cells_MOp_batches.csv",
                  label_file = "10x_cells_MOp_labels.csv")

# Train the model on the data
vae = VAE(data.nb_genes, n_batch=batches)
trainer = UnsupervisedTrainer(vae, data, train_size=0.75, use_cuda=use_cuda, frequency=5)
trainer.train(n_epochs=n_epochs, lr=lr)

ll_train = trainer.history["ll_train_set"]
ll_test = trainer.history["ll_test_set"]
x = np.linspace(1, len(ll_train), len(ll_train))
plt.plot(x, ll_train)
plt.plot(x, ll_test)
plt.ylim(min(ll_train)-50, max(ll_test) + 50)

trainer.train_set.show_t_sne(color_by='batches and labels')

latent = trainer.get_all_latent_and_imputed_values(save_latent=True,
         filename_latent=os.path.join(out, "rdim.csv"))["latent"]

imputation_dic = trainer.get_all_latent_and_imputed_values(save_imputed=True,
                 filename_imputation=os.path.join(out, "inputed.csv"))`
