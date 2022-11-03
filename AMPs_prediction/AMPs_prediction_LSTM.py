## usage python3 prediction_10lstm.py bact.txt 10lstm_bact.txt
from keras.models import load_model
from numpy import loadtxt, savetxt
from sys import argv

model = load_model('lstm.h5')
x = loadtxt('amp.txt', delimiter=",")

preds = model.predict(x)
savetxt('amp_lstm.txt', preds, fmt="%.8f", delimiter=",")