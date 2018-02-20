context("Test Functions")

test_that("align_to_unip", {
    ens.seq1 <- AAString("MDLSALREEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLK")
    ens.seq2 <- AAString("MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLA")
    ens.seq3 <- AAString("MDEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLA")
    uni.seq1 <- AAString("MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLA")
    ensembl_mut_position <- 25
    
    expect_that(align_to_unip(ens.seq1, uni.seq1, ensembl_mut_position)[[2]],
                equals(26))
    expect_that(align_to_unip(ens.seq2, uni.seq1, ensembl_mut_position)[[2]],
                equals(25))
    expect_that(align_to_unip(ens.seq3, uni.seq1, ensembl_mut_position)[[2]],
                equals(31))
    
})


test_that("calculate_boundary", {
  expect_that(calculate_boundary(1000, 10000, c(200, 400)),
              equals(list(region_1 = c(900, 1100),
                          region_2 = c(800, 1200))))
  expect_that(calculate_boundary(10, 10000, c(200, 400)),
              equals(list(region_1 = c(1, 201),
                          region_2 = c(1, 401))))
  expect_that(calculate_boundary(950, 1000, c(200, 400)),
              equals(list(region_1 = c(800, 1000),
                          region_2 = c(600, 1000))))
})

test_that("identify_hotspots",{
  
  data("SnpData", package = "DominoEffect")
  data("TestData", package = "DominoEffect")
  data("DominoData", package = "DominoEffect")
  hotspot_mutations <- identify_hotspots(mutation_dataset = TestData, 
                                         gene_data = DominoData, 
                                         snp_data = SnpData)
  hotspot_mutations <- identify_hotspots(mutation_dataset = TestData, 
                                         gene_data = DominoData, 
                                         snp_data = SnpData,
                                         flanking_region = 300)
  
  expect_that(unique(hotspot_mutations$Gene), equals(c("GNAQ", "GNA11", 
                                                       "SF3B1", "PRMT8", "CHEK2")))
  expect_that(min(hotspot_mutations$Adj_p_value_region_1), equals(1.58691e-78))
  expect_that(mean(hotspot_mutations$Adj_p_value_region_1), equals(2.006675e-13))
  expect_that(sum(hotspot_mutations$N_mut), equals(95))
  
})

test_that("DominoEffect",{
    
 
})

test_that("map_to_func_elem",{
    
    
})

