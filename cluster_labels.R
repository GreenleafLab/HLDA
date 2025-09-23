# Labels for clusters

BroadClust <- list(
    "Ep" = "Epithelial",
    "En" = "Endothelial",
    "As" = "Airway stromal",
    "Vs" = "Vascular stromal",
    "Ly" = "Lymphoid",
    "My" = "Myeloid/B-cells",
    "Le" = "Lymphatic endothelial",
    "Ch" = "Chondrocytes"
)

BroadNamedClust <- list(
    "Ep03" = "APr",
    
    "As01" = "aSMC",
    "As02" = "AirF",
    "As03" = "MyoF",
    "As04" = "AirF",
    "As05" = "AirF",
    "As06" = "AdvF",
    "Ch01" = "Chdr",
    "En01" = "gCap",
    "En02" = "Artr",
    "En03" = "Aero",
    "En04" = "Veno",
    "Ep01" = "AT2l",
    "Ep02" = "AT1l",
    
    "Ep04" = "AlvPr",
    "Ep05" = "AT1l",
    "Ep06" = "APr",
    "Ep07" = "AT1l",
    "Ep08" = "EpiC",
    "Ep09" = "Cili",
    "Ep10" = "PNEC",
    "Ep11" = "AlvPr",
    "Le01" = "Lymp",
    "Ly01" = "NK",
    "Ly02" = "Tc",
    "My01" = "APCs",
    "My02" = "Bc",
    "Vs01" = "Peri",
    "Vs02" = "vSMC"
    )

FineNamedClust <- list(
    # Epithelial
    "Ep03" = "APr",
    "Ep06" = "APr",
    "Ep09" = "Cili",
    "PNEC_2" = "PNEC1", # ASCL1 + NEUROD1+
    "PNEC_1" = "PNEC2", # NEUROD1 hi
    "PNEC_0" = "PNEC3", # ASCL1 hi. More skewed to later time points
    "Ep04" = "AlvPr",
    "Ep11" = "AlvPr",
    "Ep02" = "AT1l",
    "Ep05" = "AT1l",
    "Ep07" = "AT1l",
    "Ep01" = "AT2l",
    "Ep08" = "cycEpi",
    
    # Endothelial
    "gCap_0" = "gCap",
    "gCap_1" = "eAero",
    "En03" = "Aero",
    "Artr_0" = "Artr1",
    "Artr_1" = "Artr2",
    "Artr_2" = "ArtrTip", #ESM1 and PRND positive
    "En04" = "Veno",
    "Le01" = "Lymp",
    "gCap_2" = "cycgCap",
    
    # Stromal
    "AirF_0" = "AirF1",
    "AirF_1" = "AirF2",
    "AirF_2" = "AirF3",
    "AirF_3" = "cycAirF", #TOP2A, MKI67, etc.
    "AirF_6" = "eAdvF", # Both AirF and AdvF. PAMR1+, TWIST2-, COL14A1+, NAV3+, BRINP3+ but BMP5- and BMPER-
    "As06" = "AdvF",
    "MyoF_0" = "MyoF1", #PDGFRA+, DACH2- EYA4+ SCN7A+
    "MyoF_1" = "MyoF2", #PDGFRA+, DACH2+ EYA4+, MYH11+
    "As01" = "aSMC",
    "Ch01" = "Chdr",
    "Peri_0" = "ePeri", # SULT1E1 +
    "Peri_1" = "lPeri", #LRRTM4 + 
    "Peri_3" = "cycPeri",
    "vSMC_0" = "vSMC1", # PLN+ MYOM2+ RCAN2+
    "vSMC_1" = "vSMC2", # ITGA11+ GPC3+ more like the AirF/AdvF population

    # Immune
    "APCs_2" = "Neut1",
    "APCs_4" = "Neut2",
    "APCs_5" = "Neut3",
    "APCs_0" = "Mono1",
    "APCs_6" = "Mono2",
    "APCs_7" = "IM",
    "APCs_1" = "DC1",
    "APCs_3" = "DC2",
    "Ly01" = "NK",
    "Tc_0" = "Tc1", #SELL+
    "Tc_1" = "Tc2", # ZBTB16+
    "Bc_0" = "Plasma",
    "Bc_1" = "Bc", #MS4A1+ MME- EBF1+
    "Bc_2" = "proBc", #MS4A1- MME+ EBF1+ https://doi.org/10.1016/j.it.2022.01.003
    "APCs_8" = "MPP"
    )

CustomNamedClust <- list(
  "Epithelial" = c(
    "APr1",
    "APr2",
    "PNEC1",
    "PNEC2",
    "PNEC3",
    "Cili",
    "eAlvPr",
    "AlvPr",
    "cycEpi",
    "AT1l",
    "AT2l"
  ),
  
  "Endothelial" = c(
    "egCap",
    "cycgCap",
    "gCap",
    "eAero",
    "Aero",
    "ArtrTip",
    "Artr1",
    "Artr2",
    "Veno",
    "Lymp"
  ),
  
  "Stromal" = c(
    "AirF1",
    "AirF2",
    "AirF3",
    "cycAirF",
    "eAdvF",
    "AdvF",
    "MyoF1",
    "MyoF2",
    "aSMC",
    "ePeri",
    "cycPeri",
    "lPeri",
    "vSMC1",
    "vSMC2",
    "Chdr"
  ),
  
  "Immune" = c(
    "Tc1",
    "Tc2",
    "NK",
    "Bc",
    "Plasma",
    
    "Mono",
    "DC1",
    "DC2",
    "IM",
    "Neut1",
    "Neut2",
    
    "MPP"
  )
)

CustomNamedClust_ClusterGroups <- list(
  "Epithelial" = c(
    "APr1",
    "APr2",
    "Cili",
    "PNEC1",
    "PNEC2",
    "PNEC3",
    "eAlvPr",
    "AlvPr",
    "AT1l",
    "AT2l",
    "cycEpi"
  ),
  
  "Endothelial" = c(
    "egCap",
    "cycgCap",
    "gCap",
    "eAero",
    "Aero",
    "ArtrTip",
    "Artr1",
    "Artr2",
    "Veno",
    "Lymp"
  ),
  
  "Stromal_a" = c(
    "AirF1",
    "AirF2",
    "AirF3",
    "cycAirF",
    "eAdvF",
    "AdvF",
    "MyoF1",
    "MyoF2",
    "aSMC",
    "Chdr"
    ),
  
  "Stromal_v" = c(
    "ePeri",
    "lPeri",
    "cycPeri",
    "vSMC1",
    "vSMC2"
  ),
  
  "Immune_m" = c(
    "Mono",
    "DC1",
    "DC2",
    "IM",
    "Neut1",
    "Neut2",
    "Bc",
    "Plasma",
    "MPP"
  ),
  "Immune_l" = c(
    "NK",
    "Tc1",
    "Tc2"
  )
)

CustomNamedClust_ClusterGroups_Tree <- list(
  "Epithelial_air" = c(
    "APr1",
    "APr2",
    "Cili",
    "PNEC1",
    "PNEC2",
    "PNEC3"
  ),
  
  "Epithelial_alv" = c(
    "eAlvPr",
    "AlvPr",
    "AT1l",
    "AT2l",
    "cycEpi"
  ),
  
  "Endothelial" = c(
    "egCap",
    "cycgCap",
    "gCap",
    "eAero",
    "Aero",
    "ArtrTip",
    "Artr1",
    "Artr2",
    "Veno",
    "Lymp"
  ),
  
  "Stromal_a" = c(
    "AirF1",
    "AirF2",
    "AirF3",
    "cycAirF",
    "eAdvF",
    "AdvF",
    "MyoF1",
    "MyoF2",
    "aSMC",
    "Chdr"
  ),
  
  "Stromal_v" = c(
    "ePeri",
    "cycPeri",
    "lPeri",
    "vSMC1",
    "vSMC2"
  ),
  
  "Immune_m" = c(
    "Mono",
    "DC1",
    "DC2",
    "IM",
    "Neut1",
    "Neut2",
    "Bc",
    "Plasma",
    "MPP"
  ),
  "Immune_l" = c(
    "NK",
    "Tc1",
    "Tc2"
  )
)


fineOrder <- list(
  "Epithelial" = c("APr1", "APr2", "APr3", "APr4", "PNEC1", "PNEC2", "PNEC3", "eCili", "lCili", 
                   "TiP1", "TiP2", "AT2l", "AT1l", "EpiC"),
  "Endothelial" = c("egCap", "lgCap", "eAero", "lAero", "eArtr", "lArtr", "Veno", "Lymp"),
  "Stromal" = c("eAlvF", "lAlvF", "AdvF", "MyoF", "aSMC", "ePeri", "lPeri", "vSMC1", "vSMC2", "Chdr", "Schw", "Meso"),
  "Immune" = c("Mono1", "Mono2", "Dc", "IM", "Neut", "NK", "Tc", "Bc", "Plasma", "MPP")
)

featureSets <- list(
  "Epithelial" = c("EPCAM", "CDH1"), # Epithelial
  "Airway" = c("SOX2", "SCGB3A2"), # Airway
  "PNEC" = c("P3H2", "NRXN1", "ASCL1", "NEUROD1", "GHRL", "GRP"), #PNEC
  "Ciliated" = c("FOXJ1", "DNAH12"), #Ciliated
  "Alveolar" = c("SFTPB", "SOX9", "ETV5", "AGER", "SFTPC"), # Alveolar
  "Cycling" = c("TOP2A", "MKI67"),
  
  "Endothelial" = c("PECAM1", "CLDN5"),
  "gCaps" = c("CA4"),
  "Aerocytes" = c("EDNRB", "S100A3"),
  "Arteries" = c("GJA4", "DKK2", "FBLN5"),
  "Venous" = c("ACKR1", "PLVAP"),
  "Lymphatic" = c("CCL21", "PROX1"),
  
  "Structural" = c("COL1A1", "COL1A2"),
  "Airway" = c("PDGFRA"),
  "AlvF" = c("MEOX2", "WNT2"),
  "AdvF" = c( "PAMR1", "SERPINF1"),
  "MyoF" = c("DACH2", "EYA4"),
  "Contractile" = c("ACTA2", "TAGLN"),
  "AiwaySM" = c("HHIP", "HPSE2"),
  "Vascular" = c("PDGFRB"),
  "Peri" = c("TRPC6", "SULT1E1", "LRRTM4", "RBFOX1", "PAG1"),
  "VascSM" = c("SLIT3", "ELN",  "CSMD1", "PLN", "ITGA11"),
  "Chondrocytes" = c("ACAN", "COL2A1"),
  "Schwann" = c("CDH19", "MPZ"),
  "Mesothelial" = c("C3", "MSLN"),
  
  "Immune" = c("PTPRC"),
  "Mono1" = c("AQP9", "VCAN", "LYZ"),
  "Dc" = c("HLA-DRA", "HLA-DQA1"),
  "IM" = c("MRC1", "MSR1", "C1QB"),
  "Neut" = c("MPO", "LTF"),
  "NK" = c("NKG7", "CCL4"),
  "T_cells" = c("IL7R", "CAMK4"),
  "B_cells" = c("EBF1", "IGHM", "MS4A1"),
  "Plasma" = c("IRF4", "JCHAIN"),
  "MPP" = c("CD34", "GATA2")
)