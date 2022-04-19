import json
#%%
'TEST'
import matplotlib.pyplot as plt
import json

p = 30

name = 'tsne_'+str(p)+'.json'
with open('/Users/victorparedes/Documents/GitHub/Auto-unispring/Python/'+name, 'r') as dataFile:
    rawData = dataFile.read()
tsne_out = json.loads(rawData)
plt.scatter(*zip(*tsne_out))
plt.show()