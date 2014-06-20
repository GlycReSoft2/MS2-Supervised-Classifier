# MS2-Supervised-Classifier

Supervised Learning Classifier of MS/MS Glycopeptides


## Training and Evaluation
Each dataset (containing a single gold standard labeled dataset and two technical replicates) is used to train a Random Forest classifier. The classifier is tested on three tasks:

* Ability to label technical replicates that in turn can train a model to predicate the labels of the expert labeled dataset
* Ability to correctly predict labels of other samples' gold standard labeled dataset
