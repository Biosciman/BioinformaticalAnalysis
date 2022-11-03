## usage python3 prediction_10att.py bact.txt 10att_bact.txt

from keras.models import load_model
from numpy import loadtxt, savetxt
from Attention import Attention_layer
from sys import argv


model = load_model('att.h5', custom_objects={'Attention_layer': Attention_layer})
x = loadtxt('amp.txt', delimiter=",")

preds = model.predict(x)
savetxt('amp_att.txt', preds, fmt="%.8f", delimiter=",")
