import pickle
with open('colon_output.pkl', 'rb') as fp:
    data = pickle.load(fp)
    print data 