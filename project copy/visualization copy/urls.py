from django.urls import path

from . import views

app_name = 'visualization'
urlpatterns = [
    path('index/', views.index, name='run'),
    path('info/', views.info, name='info'),
    path('visual/', views.visual, name='visual'),
    path('', views.upload, name='upload'),
]