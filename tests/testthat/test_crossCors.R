test_that("crossCors output matches cor.test output", {
	drugActData <- exprs(getAct(rcellminerData::drugData))
  crossCors.out <- crossCors(drugActData[c("94600"), , drop=FALSE], drugActData[c("727625", "670655"), ])
  
  expect_equivalent(cor.test(as.numeric(drugActData["94600", ]), as.numeric(drugActData["727625", ]))$estimate,
                    crossCors.out$cor["94600", "727625"])
  
  expect_equivalent(cor.test(as.numeric(drugActData["94600", ]), as.numeric(drugActData["670655", ]))$estimate,
                    crossCors.out$cor["94600", "670655"])
  
  expect_equivalent(cor.test(as.numeric(drugActData["94600", ]), as.numeric(drugActData["727625", ]))$p.value,
                    crossCors.out$pval["94600", "727625"])
  
  expect_equivalent(cor.test(as.numeric(drugActData["94600", ]), as.numeric(drugActData["670655", ]))$p.value,
                    crossCors.out$pval["94600", "670655"])
})

test_that("crossCorsSpearman computations are correct", {
	nX <- 4
	nY <- 5
	
	drugActData <- exprs(getAct(rcellminerData::drugData))
	X <- drugActData[1:nX, ]
	Y <- drugActData[1000:(1000 + nY - 1), ]
	
	xcorResults <- crossCors(X, Y, method = "spearman")
	
	for (xName in rownames(xcorResults$cor)){
		for (yName in colnames(xcorResults$cor)){
			xvec <- as.numeric(drugActData[xName, ])
			yvec <- as.numeric(drugActData[yName, ])
			
			expectedCor <- cor(xvec, yvec, use = "pairwise.complete.obs", method = "spearman")
			expect_identical(expectedCor, xcorResults$cor[xName, yName])
			
			expectedPval <- suppressWarnings(cor.test(xvec, yvec, method = "spearman")$p.value)
			expect_identical(expectedPval, xcorResults$pval[xName, yName])
		}
	}
	
	xcorResults <- crossCors(X, method = "spearman")
	
	for (xName in rownames(xcorResults$cor)){
		for (yName in colnames(xcorResults$cor)){
			xvec <- as.numeric(drugActData[xName, ])
			yvec <- as.numeric(drugActData[yName, ])
			
			expectedCor <- cor(xvec, yvec, use = "pairwise.complete.obs", method = "spearman")
			expect_identical(expectedCor, xcorResults$cor[xName, yName])
			
			expectedPval <- suppressWarnings(cor.test(xvec, yvec, method = "spearman")$p.value)
			expect_identical(expectedPval, xcorResults$pval[xName, yName])
		}
	}
	
})