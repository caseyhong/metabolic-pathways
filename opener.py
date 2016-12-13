import pickle
with open('lung_output.pkl', 'rb') as fp:
    data = pickle.load(fp)
    print data 