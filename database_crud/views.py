from rest_framework import viewsets
from rest_framework.views import APIView
from rest_framework.decorators import action
from .pagination import MoleculePagination
from django.http import HttpResponse
from .models import Molecule
import json
import csv
from .serializers import MoleculeSerializer
from rest_framework.response import Response
from rest_framework.parsers import MultiPartParser
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors, Crippen, Lipinski


class MoleculeViewSet(viewsets.ModelViewSet):
    queryset = Molecule.objects.all()
    serializer_class = MoleculeSerializer

    @action(detail=False, methods=['get'])
    def lipinski_filter(self, request):
        """
        Filters molecules based on Lipinski's Rule of Five.
        Allows selecting molecules that match at least N criteria.
        """
        min_criteria = int(request.query_params.get(
            'min_criteria', 4))  # Default: At least 3 criteria

        # Get ordering field from query params (default is 'id')
        order_by = request.query_params.get('order_by', 'id')

        # Get order direction (default is ascending)
        order = request.query_params.get('order', 'asc')

        allowed_fields = [field.name for field in Molecule._meta.fields]

        if order_by not in allowed_fields:
            return Response({"error": f"Invalid field. Allowed fields: {allowed_fields}"}, status=400)

        # Get all molecules
        molecules = Molecule.objects.all()
        filtered_molecules = []

        for molecule in molecules:
            criteria_met = 0
            if molecule.mol_weight <= 500:
                criteria_met += 1
            if molecule.clogp <= 5:
                criteria_met += 1
            if molecule.hb_donors <= 5:
                criteria_met += 1
            if molecule.hb_acceptors <= 10:
                criteria_met += 1

            if criteria_met >= min_criteria:
                filtered_molecules.append(molecule)

        if order == 'asc':
            filtered_molecules.sort(key=lambda x: getattr(x, order_by))
        else:
            filtered_molecules.sort(
                key=lambda x: getattr(x, order_by), reverse=True)

        page = self.paginate_queryset(filtered_molecules)
        if page is not None:
            serializer = self.get_serializer(page, many=True)
            return self.get_paginated_response(serializer.data)

        # Serialize results
        serializer = self.get_serializer(filtered_molecules, many=True)
        return Response(serializer.data)

    @action(detail=False, methods=['get'])
    def lead_like_filter(self, request):
        """
        Filters molecules based on Lead-like Rule of three.
        Allows selecting molecules that match at least N criteria.
        """
        min_criteria = int(request.query_params.get(
            'min_criteria', 5))  # Default: At least 3 criteria

        # Get ordering field from query params (default is 'id')
        order_by = request.query_params.get('order_by', 'id')
        # Get order direction (default is ascending)
        order = request.query_params.get('order', 'asc')
        allowed_fields = [field.name for field in Molecule._meta.fields]

        if order_by not in allowed_fields:
            return Response({"error": f"Invalid field. Allowed fields: {allowed_fields}"}, status=400)

        # Get all molecules
        molecules = Molecule.objects.all()
        filtered_molecules = []

        for molecule in molecules:
            criteria_met = 0
            if molecule.mol_weight <= 300:
                criteria_met += 1
            if molecule.clogp <= 3:
                criteria_met += 1
            if molecule.hb_donors <= 3:
                criteria_met += 1
            if molecule.hb_acceptors <= 3:
                criteria_met += 1
            if molecule.num_rot_bonds <= 3:
                criteria_met += 1

            if criteria_met >= min_criteria:
                filtered_molecules.append(molecule)
        if order == 'asc':
            filtered_molecules.sort(key=lambda x: getattr(x, order_by))
        else:
            filtered_molecules.sort(
                key=lambda x: getattr(x, order_by), reverse=True)

        page = self.paginate_queryset(filtered_molecules)
        if page is not None:
            serializer = self.get_serializer(page, many=True)
            return self.get_paginated_response(serializer.data)

        # Serialize results
        serializer = self.get_serializer(filtered_molecules, many=True)
        return Response(serializer.data)

    @action(detail=False, methods=['get'])
    def bioaviability_filter(self, request):
        """
        Filters molecules based on bioaviability.
        Allows selecting molecules that match at least N criteria.
        """
        min_criteria = int(request.query_params.get(
            'min_criteria', 6))  # Default: At least 3 criteria

        order_by = request.query_params.get('order_by', 'id')

        # Get order direction (default is ascending)
        order = request.query_params.get('order', 'asc')

        allowed_fields = [field.name for field in Molecule._meta.fields]
        print(order_by)
        if order_by not in allowed_fields:
            return Response({"error": f"Invalid field. Allowed fields: {allowed_fields}"}, status=400)

        # Get all molecules
        molecules = Molecule.objects.all()
        filtered_molecules = []

        for molecule in molecules:
            criteria_met = 0
            if molecule.mol_weight <= 500:
                criteria_met += 1
            if molecule.clogp <= 5:
                criteria_met += 1
            if molecule.hb_donors <= 5:
                criteria_met += 1
            if molecule.hb_acceptors <= 10:
                criteria_met += 1
            if molecule.tpsa <= 200:
                criteria_met += 1
            if molecule.num_rings_alt <= 10:
                criteria_met += 1
            if molecule.num_rot_bonds <= 5:
                criteria_met += 1

            if criteria_met >= min_criteria:
                filtered_molecules.append(molecule)
        if order == 'asc':
            filtered_molecules.sort(key=lambda x: getattr(x, order_by))
        else:
            filtered_molecules.sort(
                key=lambda x: getattr(x, order_by), reverse=True)

        # Apply pagination (if defined)
        page = self.paginate_queryset(filtered_molecules)
        if page is not None:
            serializer = self.get_serializer(page, many=True)
            return self.get_paginated_response(serializer.data)

        # If pagination is not set, return the full list
        serializer = self.get_serializer(filtered_molecules, many=True)
        return Response(serializer.data)

    @action(detail=False, methods=['get'])
    def order_by(self, request):
        """
        Custom endpoint to return molecules ordered by any field.
        Use ?order_by=field_name and ?order=desc for descending order.
        """
        # Get ordering field from query params (default is 'id')
        order_by = request.query_params.get('order_by', 'id')

        # Get order direction (default is ascending)
        order = request.query_params.get('order', 'asc')

        # List of allowed fields to prevent SQL injection
        allowed_fields = [field.name for field in Molecule._meta.fields]

        if order_by not in allowed_fields:
            return Response({"error": f"Invalid field. Allowed fields: {allowed_fields}"}, status=400)

        # Determine sorting order
        ordering = f"-{order_by}" if order == "desc" else order_by

        # Query the database
        molecules = Molecule.objects.order_by(ordering)

        # Apply pagination (if defined)
        page = self.paginate_queryset(molecules)
        if page is not None:
            serializer = self.get_serializer(page, many=True)
            return self.get_paginated_response(serializer.data)

        # If pagination is not set, return the full list
        serializer = self.get_serializer(molecules, many=True)
        return Response(serializer.data)

    @action(detail=False, methods=['get'])
    def get_atoms_coordinates_and_bonds(self, request):
        """
        Custom endpoint to get the atoms coordinates and bonds of a molecule.
        requires the molecule id as a parameter.
        """
        molecule_id = request.query_params.get('id')

        if not molecule_id:
            return Response({"error": "Missing molecule_id parameter"}, status=400)

        try:
            molecule = Molecule.objects.get(id=molecule_id)
        except Molecule.DoesNotExist:
            return Response({"error": "Molecule not found"}, status=404)

        # Get the molecule's SMILES string
        smiles = molecule.smiles
        mol = Chem.MolFromSmiles(smiles)
        AllChem.Compute2DCoords(mol)
        # Get the atoms coordinates
        conf = mol.GetConformer()

        # Extract atomic coordinates
        atoms = {atom.GetIdx(): {
            "symbol": atom.GetSymbol(),
            "x": conf.GetAtomPosition(atom.GetIdx()).x,
            "y": conf.GetAtomPosition(atom.GetIdx()).y
        } for atom in mol.GetAtoms()}

        # Extract bonds
        bonds = [{"start": bond.GetBeginAtomIdx(),
                  "end": bond.GetEndAtomIdx(),
                  "type": str(bond.GetBondType())} for bond in mol.GetBonds()]

        return Response({"atoms": atoms, "bonds": bonds})

    @action(detail=False, methods=['get'])
    def search_by_name_or_formula(self, request):
        """
        Custom endpoint to search molecules by name or formula.
        """
        query = request.query_params.get('query')

        if not query:
            return Response({"error": "Missing query parameter"}, status=400)

        # Search by name or formula
        molecules = Molecule.objects.filter(
            name__icontains=query) | Molecule.objects.filter(formula__icontains=query)

        # Serialize results
        serializer = self.get_serializer(molecules, many=True)
        return Response(serializer.data)

class SDFUploadView(APIView):
    parser_classes = [MultiPartParser]

    def post(self, request, *args, **kwargs):
        sdf_file = request.FILES.get('file')

        if not sdf_file:
            return Response({"error": "No file provided"}, status=400)

        # Read the SDF file content
        sdf_content = sdf_file.read().decode('utf-8')

        suppl = Chem.SDMolSupplier()
        suppl.SetData(sdf_content)

        molecules = []
        for mol in suppl:
            if mol is None:
                continue  # Skip invalid molecules

            try:
                molecular_id = mol.GetProp('Mol_ID') if mol.HasProp(
                    "Mol_ID") else "Unknown"
                name = mol.GetProp("_Name") if mol.HasProp(
                    "_Name") else "Unknown"
                smiles = Chem.MolToSmiles(mol)
                hb_donors = Lipinski.NumHDonors(mol)
                hb_acceptors = Lipinski.NumHAcceptors(mol)
                clogp = Crippen.MolLogP(mol)
                num_rot_bonds = Lipinski.NumRotatableBonds(mol)
                num_rings_alt = Chem.rdMolDescriptors.CalcNumRings(mol)
                tpsa = Descriptors.TPSA(mol)
                mol_weight = Descriptors.MolWt(mol)
                formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                # Approximate LogD calculation
                logD = float(mol.GetProp('LogD'))

                molecule = Molecule(
                    molecular_id=molecular_id,
                    name=name,
                    smiles=smiles,
                    hb_donors=hb_donors,
                    hb_acceptors=hb_acceptors,
                    clogp=clogp,
                    num_rot_bonds=num_rot_bonds,
                    num_rings_alt=num_rings_alt,
                    tpsa=tpsa,
                    mol_weight=mol_weight,
                    formula=formula,
                    logD=logD,
                )
                molecules.append(molecule)

            except Exception as e:
                print(f"Error processing molecule: {e}")

        # Bulk insert into database
        Molecule.objects.bulk_create(molecules)

        return Response({"message": f"{len(molecules)} molecules uploaded successfully!"}, status=201)
