module load python/3.6.4_gcc5_np1.14.5
module load cuda/8.0
module load pytorch/0.1.5

pip install -t ~/scVI/ --user --upgrade-strategy="only-if-needed" 
