N1 = 2e5
N2_seq = c(2e4,2e5)
N2_label = c("20k", "200k")

scenarios = c("1" = "5 shared causal variants", 
              "2" = "3 shared causal variants",
              "3" = "1 shared causal variant")
lds = c("In-sample LD", "1kG LD")

methods = c("MACHINE + g-LDSC", "MACHINE + PolyFun", "MACHINE",
            "MESuSiE", "SuSiEx", "XMAP", "MultiSuSiE", 
            "h2-D2 + g-LDSC", "h2-D2 + PolyFun", "h2-D2",
            "SuSiE", "RSparsePro", "CARMA")

ld_methods = list(
  "In-sample LD" = c(
    "MACHINE + g-LDSC" = "MACHINE-gLDSC",
    "MACHINE + PolyFun" = "MACHINE-polyfun",
    "MACHINE" = "MACHINE",
    "MESuSiE" = "MESuSiE",
    "SuSiEx" = "SuSiEx",
    "XMAP" = "XMAP",
    "MultiSuSiE" = "MultiSuSiE",
    "h2-D2 + g-LDSC" = "h2D2-gLDSC",
    "h2-D2 + PolyFun" = "h2D2-polyfun",
    "h2-D2" = "h2D2",
    "SuSiE" = "SuSiE",
    "CARMA" = "CARMA"),
  "1kG LD" = c(
    "MACHINE + g-LDSC" = "MACHINE_1kg-gLDSC",
    "MACHINE + PolyFun" = "MACHINE_1kg-polyfun",
    "MACHINE" = "MACHINE_1kg",
    "h2-D2 + g-LDSC" = "h2D2_1kg-gLDSC",
    "h2-D2 + PolyFun" = "h2D2_1kg-polyfun",
    "h2-D2" = "h2D2_1kg",
    "RSparsePro" = "RSparsePro_1kg",
    "CARMA" = "CARMA_1kg"))

pops = c("cross" = "Cross", "pop1" = "EUR", "pop2" = "EAS", "shared" = "Shared")

n_causal = data.frame("Cross" = c(5,5,5), 
                      "EUR" = c(5,4,3), 
                      "EAS" = c(5,4,3),
                      "Shared" = c(5,3,1))