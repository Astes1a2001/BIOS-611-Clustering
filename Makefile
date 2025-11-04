RSCRIPT = Rscript

# default target
all: figs/kmeans_hypercube_gap.png figs/spectral_shells_gap.png

# ensure folders exist
results:
	mkdir -p results

figs:
	mkdir -p figs

# K-means on hypercube data
figs/kmeans_hypercube_gap.png: scripts/01_kmeans_hypercube.R figs results
	$(RSCRIPT) scripts/01_kmeans_hypercube.R

# Spectral clustering on shells
figs/spectral_shells_gap.png: scripts/02_spectral_helpers.R scripts/03_spectral_shells.R figs results
	$(RSCRIPT) scripts/03_spectral_shells.R

# Remove generated outputs
clean:
	rm -rf figs results
