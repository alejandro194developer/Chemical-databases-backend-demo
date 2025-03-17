from rest_framework.pagination import PageNumberPagination

class MoleculePagination(PageNumberPagination):
    page_size = 5  # Number of molecules per page
    page_size_query_param = 'page_size'  # Allow clients to set page size
    max_page_size = 50  # Max items per page to prevent overload
