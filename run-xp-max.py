#!/usr/bin/env python3

from graphTools import *
from expTools import *
import os

# Dictionnaire avec les options de compilations d'apres commande
options = {}
options["-s "] = [2048]
options["-k "] = ["max"]
options["-v "] = ["omp","task","ordered"]
options["-a "] = [50,125,250]
options["-l "] = ["images/spirale.png"]
options["-nt "] = [16]

# Pour renseigner l'option '-of' il faut donner le chemin depuis le fichier easypap
#options["-of "] = ["./plots/data/perf_data.csv"]


# Dictionnaire avec les options OMP
ompenv = {}
#ompenv["OMP_NUM_THREADS="] = [1] + list(range(2, 9, 2))
#ompenv["OMP_PLACES="] = ["cores", "threads"]
#ompenv["OMP_SCHEDULE="] = ["dynamic","static","static,1"]
nbrun = 10
# Lancement des experiences
execute('./run ', ompenv, options, nbrun, verbose=True, easyPath=".")

# Lancement de la version seq avec le nombre de thread impose a 1
options["-v "] = ["seq"]
ompenv["OMP_NUM_THREADS="] = [1]
execute('./run', ompenv, options, nbrun, verbose=False, easyPath=".")