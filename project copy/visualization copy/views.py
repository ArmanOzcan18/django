from django.shortcuts import render
from django.http import HttpResponse
from django.http import HttpResponseRedirect
from django.urls import reverse
from .forms import UploadForm
from .algorithms import plot
from .algorithms import triangles
from visualization.data import Data
import logging
import pickle

def index(request):
    if request.method == 'POST' and 'run_script' in request.POST:
        global validated_triangles, removed_triangles, removed_seacells
        global adjacency_matrix, celltypes, UMAP_coords
        # import function to run
        logging.debug(filename)

        with open(filename, 'rb') as f:
            input = pickle.load(f)  # deserialize using load()
        f.close()

        data = triangles(input.model, input.ad)
        validated_triangles = data.validated_triangles
        removed_triangles = data.removed_triangles
        removed_seacells = data.removed_seacells
        adjacency_matrix = data.adjacency_matrix
        celltypes = data.celltypes
        UMAP_coords = data.UMAP_coords
        # return user to required page
        return HttpResponseRedirect(reverse("visualization:info"))
    else:
        return render(request, "index.html")

def info(request):
    if request.method == 'POST' and 'visualize' in request.POST:
        # return user to required page
        return HttpResponseRedirect(reverse("visualization:visual"))
    else:
        context = {'validated_triangles': validated_triangles, 'removed_triangles': removed_triangles, 'len': len(removed_seacells), }
        return render(request, "info.html", context)

def visual(request):
    html = plot(adjacency_matrix, celltypes, UMAP_coords)
    return HttpResponse(html)
def upload(request):
    if request.method == 'POST':
        form = UploadForm(request.POST, request.FILES)
        if form.is_valid():
            # Save the uploaded file to the database
            upload = form.save()
            # Get the file name
            global filename
            filename = upload.file.name.split('/')[-2] + "/" + upload.file.name.split('/')[-1]
            sent = {
                'file_name': upload.file.name.split('/')[-1],
            }
            return HttpResponseRedirect(reverse('visualization:run'))
    else:
        form = UploadForm()
    return render(request, 'upload.html', {'form': form})