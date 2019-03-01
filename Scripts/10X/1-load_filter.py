from scvi.dataset import LoomDataset, CsvDataset
import urllib.request
import os
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from scvi.models import *
from scvi.inference import UnsupervisedTrainer

# Set the parameters
loc="/pylon5/ib5phhp/hectorrb/10x_cells_MOp/"
out="/pylon5/ib5phhp/hectorrb/output/
n_epochs_all = None
show_plot = True
n_epochs=400 if n_epochs_all is None else n_epochs_all
lr=1e-3
use_batches=False
use_cuda=True

# Load the data
data = CsvDataset("_umi.csv.gz",
                  save_path = loc,
                  compression='gzip')

# Train the data
vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * use_batches)
trainer = UnsupervisedTrainer(vae, gene_dataset, train_size=0.75,
                              use_cuda=use_cuda, frequency=5)
trainer.train(n_epochs=n_epochs, lr=lr)

ll_train = trainer.history["ll_train_set"]
ll_test = trainer.history["ll_test_set"]
x = np.linspace(0,50,(len(ll_train)))
plt.plot(x, ll_train)
plt.plot(x, ll_test)
plt.ylim(min(ll_train)-50, max(ll_test) + 50)

trainer.train_set.show_t_sne(color_by='batches and labels')

latent = trainer.get_all_latent_and_imputed_values(save_latent=True,
         filename_latent=os.path.join(out, "rdim.csv"))["latent"]

imputation_dic = trainer.get_all_latent_and_imputed_values(save_imputed=True,
                 filename_imputation=os.path.join(out, "inputed.csv"))`
