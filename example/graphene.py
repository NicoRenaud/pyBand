import pyBand


bs = pyBand.BandStructure('graphene.in','parameters')
bs.electronic_structure()
bs.compute_bands()
bs.pickle()
bs.plot_bands()