from django.urls import path, include
from rest_framework.routers import DefaultRouter
from .views import MoleculeViewSet
from .views import SDFUploadView

router = DefaultRouter()
router.register(r'molecules', MoleculeViewSet)

urlpatterns = [
    path('', include(router.urls)),
    path('data/upload-sdf/', SDFUploadView.as_view(), name='upload_sdf'),
]
