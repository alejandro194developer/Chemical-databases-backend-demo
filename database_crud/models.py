from django.db import models

class Molecule(models.Model):
    molecular_id = models.CharField(max_length=255)
    name = models.TextField()
    smiles = models.TextField()
    hb_donors = models.IntegerField()
    hb_acceptors = models.IntegerField()
    clogp = models.FloatField()
    num_rot_bonds = models.IntegerField()
    num_rings_alt = models.IntegerField()
    tpsa = models.FloatField()
    mol_weight = models.FloatField()
    formula = models.TextField()
    logD = models.FloatField()

    def __str__(self):
        return f"{self.name} ({self.molecular_id})"
