from django.shortcuts import render
from django.http import HttpResponse

import pickle

with open('/Users/armanozcan/Desktop/project/visualization/data.pkl', 'rb') as f:
    data = pickle.load(f) # deserialize using load()
f.close()

validated_triangles = data.validated_triangles
removed_triangles = data.removed_triangles
removed_seacells = data.removed_seacells

def index(request):
    html = "<html><body><br><br><h1> Welcome to the Main Page!</h1><h3>Here are the data taken from our model:</h3><ul><li>The set of confirmed triangles are: <br> %s </li><br><li>The set of removed triangles are: <br> %s. </li><br><li>The Number of removed SEACELLS is %s. </li></ul></body></html>" % (validated_triangles, removed_triangles, len(removed_seacells))
    return HttpResponse(html)

