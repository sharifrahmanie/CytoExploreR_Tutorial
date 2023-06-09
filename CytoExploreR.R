capabilities() 
# If tcltk = FALSE 
# Installing R from source for Linux machines
# https://cran.r-project.org/src/base/R-4/R-4.2.3.tar.gz
# sudo apt-get install tcl-dev tk-dev
# ./configure --enable-R-shlib=yes --with-tcltk=/usr/lib/tcltk --with-x
# make 
# sudo make install 

# For windows users
# https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html
# https://www.activestate.com/products/tcl/

capabilities() 

install.packages("remotes")
# Bioconductor
install.packages("BiocManager")
# Install cytoinstaller
remotes::install_github("RGLab/cytoinstaller")
# Install cytoverse packages

remotes::install_github("RGLab/cytoverse")
# Update some missing dependencies
cytoverse::cytoverse_update()
# CytoExploreRData 
remotes::install_github("DillonHammill/CytoExploreRData")
# CytoExploreR 
remotes::install_github("DillonHammill/CytoExploreR")
# If it doesn't work
# Download https://github.com/DillonHammill/CytoExploreRData.git
# R CMD build --no-build-vignettes CytoExploreRData


# Load required packages
require(CytoExploreR)
require(CytoExploreRData)

help(package = "CytoExploreR")
# Compensation FCS Files
cyto_save(Compensation, 
          save_as = "Compensation-Samples")

# Activation FCS Files
cyto_save(Activation,
          save_as = "Activation-Samples")

# Setup compensation controls
gs <- cyto_setup("Compensation-Samples",
                 gatingTemplate = "Compensation-gatingTemplate.csv")

# Transform fluorescent channels - default logicle transformations
gs <- cyto_transform(gs)

# Gate Cells
cyto_gate_draw(gs,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A", "SSC-A"))

# Gate Single Cells
cyto_gate_draw(gs,
               parent = "Cells",
               alias = "Single Cells",
               channels = c("FSC-A", "FSC-H"))

# Compute spillover matrix
spill <- cyto_spillover_compute(gs,
                                parent = "Single Cells",
                                spillover = "Spillover-Matrix.csv")
# Edit spillover matrix
spill <- cyto_spillover_edit(gs,
                             parent = "Single Cells",
                             spillover = "Spillover-Matrix.csv")
# Visualise uncompensated data
cyto_plot_compensation(gs,
                       parent = "Single Cells")

# Visualise compensated data
cyto_plot_compensation(gs,
                       parent = "Single Cells",
                       spillover = "Spillover-Matrix.csv",
                       compensate = TRUE)

# Analyse Samples
# Load and annotate samples
gs <- cyto_setup("Activation-Samples",
                 gatingTemplate = "Activation-gatingTemplate.csv")
# Apply compensation
gs <- cyto_compensate(gs,
                      spillover = "Spillover-Matrix.csv")

# Transform fluorescent channels - default logicle transformations
gs <- cyto_transform(gs)

# Gate Cells
cyto_gate_draw(gs,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"))
# Gate Single Cells
cyto_gate_draw(gs,
               parent = "Cells",
               alias = "Single Cells",
               channels = c("FSC-A","FSC-H"))
# Extract unstained control 
NIL <- cyto_extract(gs, "Single Cells")[[33]]

# Gate Live Cells
cyto_gate_draw(gs,
               parent = "Single Cells",
               alias = c("Dead Cells", "Live Cells"),
               channels = c("Hoechst-405", "Hoechst-430"),
               type = "rectangle",
               negate = TRUE,
               overlay = NIL)
# Gate T Cells and Dedritic Cells
cyto_gate_draw(gs,
               parent = "Live Cells",
               alias = c("T Cells", "Dendritic Cells"),
               channels = c("CD11c", "Va2"),
               type = c("ellipse", "rectangle"),
               overlay = NIL)
# Gate CD4 & CD8 T Cells
cyto_gate_draw(gs,
               parent = "T Cells",
               alias = c("CD4 T Cells", "CD8 T Cells"),
               channels = c("CD4", "CD8"),
               type = "r")
# Extract CD4 T Cells
CD4 <- cyto_extract(gs, 
                    parent = "T Cells",
                    select = "CD4 T Cells")
OVA_0 <- cyto_select(Activation,
                  OVAConc = 0)

# Extract naive CD4 T Cells
CD4_naive <- CD4[cyto_details(OVA_0)[[1]]]

# Gate CD69+ CD4 T Cells
cyto_gate_draw(gs,
               parent = "CD4 T Cells",
               alias = "CD69+ CD4 T Cells",
               channels = c("CD44", "CD69"),
               type = "rectangle",
               overlay = CD4_naive)

# Gate CD69+ CD8 T Cells
cyto_gate_draw(gs,
               parent = "CD8 T Cells",
               alias = "CD69+ CD8 T Cells",
               channels = c("CD44", "CD69"),
               type = "rectangle",
               contour_lines = 15)

# Gating Tree
cyto_plot_gating_tree(gs[[32]],
                      stat = "freq")
# Gating scheme
cyto_plot_gating_scheme(gs[32],
                        back_gate = TRUE,
                        gate_track = TRUE)
# Compute medFI - exclude unstained control
cyto_stats_compute(gs[1:32],
                   alias = c("CD69+ CD4 T Cells",
                             "CD69+ CD8 T Cells"),
                   stat = "median",
                   channels = c("CD44", "CD69"))
