#' Subfunction to automatically draw up custom palettes for BGCs.
#'
#' Takes user guidelines about R color package and module and the number of gene-types to be highlighted, and designs a palette for the genes.  Adds light/dark colors if over a dozen colors are required, and by default highlights hypothetical proteins ("hyp"), mobile genetic elements ("mge"), and genes not specified for highlighting by the user ("other") with grayscale colors.
#' @param colorType R color package from which palettes can be drawn.
#' @param paletteInput Specific color scheme within your desired R color package.
#' @param processed Partially processed data frame from IMG or analyzeNeighbors with metadata for your gene (family) of interest.  Has already been through the initial prettyClusterDiagrams input process.
#' @return List containing sub-list with custom palette and sub-list with gene-types
#' @export
#' @importFrom utils read.csv write.csv write.table read.table
#' @importFrom rlang .data
#' @importFrom magrittr %>%   
#' @examples 
#' \dontrun{
#' pcdColor <- pcdColor(colorType = colorType, 
#'                                    paletteInput = paletteInput, 
#'                                    processed = processed)
#' }
#'
 pcdColor <- function(colorType = colorType, paletteInput = paletteInput, processed = processed) {
        ## this will color genes according the the genes you specified
        ## in its default form it has distinct colors for hypothetical proteins, unspecified proteins, and mobile genetic elements
        ## it will also add a 2nd light/dark color shade if there are > 10 named genes
        ## standardizing on a minimal palette set so that all get colorRamped similarly
        ## fish!
        if (colorType == "fishualize") {
            if(requireNamespace("fishualize", quietly=TRUE)) {
                tempPalette <- fishualize::fish(5,option=paletteInput)
            } else {
                print("Chosen palette not installed, using viridis.")
                tempPalette <- viridis::viridis(5,option=paletteInput)
            }
        }
        if (colorType == "ghibli") {
            if(requireNamespace("ghibli", quietly=TRUE)) {
                tempPalette <- ghibli::ghibli_palette(paletteInput)[1:7]
            } else {
                print("Chosen palette not installed, using viridis.")
                tempPalette <- viridis::viridis(5,option=paletteInput)
            }
        }
        if (colorType == "lisa") { 
            if(requireNamespace("lisa", quietly=TRUE)) {
                tempPalette <- lisa::lisa_palette(paletteInput, 5)[1:5]
            } else {
                print("Chosen palette not installed, using viridis.")
                tempPalette <- viridis::viridis(5,option=paletteInput)
            }
        }
        if (colorType == "nord") {
            if(requireNamespace("nord", quietly=TRUE)) {
                tempPalette <- nord::nord(paletteInput,5)      
            } else {
                print("Chosen palette not installed, using viridis.")
                tempPalette <- viridis::viridis(5,option=paletteInput)
            }
        }
        if (colorType == "rtist") {
            if(requireNamespace("rtist", quietly=TRUE)) {
                tempPalette <- rtist::rtist_palette(paletteInput,5)[1:5]      
            } else {
                print("Chosen palette not installed, using viridis.")
                tempPalette <- viridis::viridis(5,option=paletteInput)
            }
        }
        if (colorType == "scico") {
            if(requireNamespace("scico", quietly=TRUE)) {
                tempPalette <- scico::scico(5,palette=paletteInput)      
            } else {
                print("Chosen palette not installed, using viridis.")
                tempPalette <- viridis::viridis(5,option=paletteInput)
            }
        }
        if (colorType == "wesanderson") {
            if(requireNamespace("wesanderson", quietly=TRUE)) {
                tempPalette <- wesanderson::wes_palette(paletteInput,5)      
            } else {
                print("Chosen palette not installed, using viridis.")
                tempPalette <- viridis::viridis(5,option=paletteInput)
            }
        }
        if (colorType == "viridis") {tempPalette <- viridis::viridis(5,option=paletteInput)}
        ## now tailoring palettes to our gene set
        geneTypes <- unique(processed$gene)
        colorNum <- length(geneTypes)
        geneTypes <- geneTypes[stringr::str_order(geneTypes)]
        allGeneTypes <- geneTypes
        notMe <- which(geneTypes =="hyp")
        notMeEither <- which(geneTypes == "mge")
        norMe <- which(geneTypes == "other")
        if (any(notMe>0, notMeEither>0, norMe>0)) {
            geneTypes <- geneTypes[-c(notMe,notMeEither,norMe)]
        }
        numGeneTypes <- length(geneTypes)
        halfGeneTypes <- numGeneTypes
        if(numGeneTypes >= 10) {
            halfGeneTypes <- ceiling(numGeneTypes/2)
        }
        if (halfGeneTypes == numGeneTypes) {
            ## keep a reasonably bright version if possible
            geneColors <- colorRampPalette(tempPalette)(numGeneTypes)
        } else if (halfGeneTypes != numGeneTypes) {
            if (colorType == "ghibli") {
                ## it has built-in medium, light, and dark options
                paletteDark <- gsub("Medium","Dark",paletteInput)
                paletteDark <- ghibli::ghibli_palette(paletteDark)
                paletteLight <- gsub("Medium","Light",paletteInput)
                paletteLight <- ghibli::ghibli_palette(paletteLight)
                darkColors <- colorRampPalette(paletteDark)(halfGeneTypes)
                lightColors <- colorRampPalette(paletteLight)(halfGeneTypes)
                geneColors <- c(darkColors, lightColors)
            } else {
                ## for the rest we roll our own, i guess
                ## these are sorta arbitrary functions
                geneColors <- colorRampPalette(tempPalette)(halfGeneTypes)
                darkColors <- rgb2hsv(col2rgb(geneColors))
                darkColors["s",] <- abs(0.7*darkColors["s",])
                darkColors["v",] <- abs(darkColors["v",] - 0.4*(darkColors["v",]))
                darkColors <- apply(darkColors,2,function(x) hsv(x[1], x[2], x[3]))
                lightColors <- rgb2hsv(col2rgb(geneColors))
                lightColors["s",] <- abs(0.4*lightColors["s",])
                lightColors["v",] <- abs(lightColors["v",] + 0.7*(1-lightColors["v",]))
                lightColors <- apply(lightColors,2,function(x) hsv(x[1], x[2], x[3]))
                geneColors <- c(darkColors, lightColors)
            }
        }       
        ## this specifically codes (unidentified) hypothetical proteins, mobile genetic elements, and known proteins not on the annotation list
        ## and reinserts them into the palette, so now it doesn't sacrifice pretty colors for them either
        ## as a default: hyp. = white, MGE = grey, 
            ## troubleshooting CASES
            ## CASE A has hyp, mge, other
                ## OK   geneTypes9 <- list("abc","hyp","mge","other","qrs")
                ## OK   geneTypes1 <- list("abc","hyp","jkl","mge","nop","other","qrs")
                ## OK   geneTypes2 <- list("abc","hyp","jkl","mge","other","qrs")
                ## OK   geneTypes8 <- list("abc","hyp","mge","nop","other","qrs")
                ## OK   geneTypes21 <- list("abc","hyp","mge","other")
                ## OK   geneTypes22 <- list("abc","hyp","jkl","mge","nop","other")
                ## OK   geneTypes23 <- list("abc","hyp","jkl","mge","other")
                ## OK   geneTypes24 <- list("abc","hyp","mge","nop","other")
                ## OK   geneTypes29 <- list("hyp","mge","other","qrs")
                ## OK   geneTypes30 <- list("hyp","jkl", "mge","other","qrs")
                ## OK   geneTypes31 <- list("hyp","jkl", "mge","nop","other","qrs")
                ## OK   geneTypes33 <- list("hyp", "mge","nop","other","qrs")
                ## OK   geneTypes32 <- list("hyp", "jkl", "mge","nop","other")
            ## CASE B has hyp, mge, no other
                ## OK   geneTypes3 <- list("abc","hyp", "jkl", "mge", "nop","qrs")
                ## OK   geneTypes13 <- list("abc","hyp","mge", "nop","qrs")
                ## OK   geneTypes19 <- list("abc","hyp", "jkl", "mge")
                ## OK   geneTypes20 <- list("abc","hyp", "mge")
                ## OK   geneTypes34 <- list("hyp", "jkl", "mge", "nop","qrs")
                ## OK   geneTypes35 <- list("hyp", "jkl", "mge")
                ## OK   geneTypes36 <- list("hyp", "mge")
            ## CASE C: has hyp, other (no mge)
                ## OK   geneTypes4 <- list("abc","hyp","jkl","nop","other","qrs")
                ## OK   geneTypes5 <- list("abc","hyp","other","qrs")
                ## OK   geneTypes17 <- list("abc","hyp","jkl","nop","other")
                ## OK   geneTypes18 <- list("abc","hyp","other")
                ## OK   geneTypes37 <- list("hyp","jkl","nop","other","qrs")
                ## OK   geneTypes38 <- list("hyp","jkl","nop","other")
                ## OK   geneTypes39 <- list("hyp","other")
            ## CASE D: hyp only
                ## OK   geneTypes25 <- list("abc","hyp", "jkl", "nop", "qrs")
                ## OK   geneTypes16 <- list("abc","hyp")
                ## OK   geneTypes40 <- list("hyp", "jkl", "nop", "qrs")
                ## OK   geneTypes41 <- list("mge","nop","other","qrs")
                ## OK geneTypes42 <- list("mge","nop","other")
                ## OK geneTypes43 <- list("mge","other")
            ## CASE E no hyp, has mge and other
                ## OK   geneTypes6 <- list("abc","jkl","mge","nop","other","qrs")
                ## OK   geneTypes7 <- list("abc","jkl","mge","other","qrs")
                ## OK   geneTypes14 <- list("abc","jkl","mge","other")
                ## OK   geneTypes15 <- list("abc","jkl","mge","nop","other")
                ## OK   geneTypes46 <- list("mge","nop","other","qrs")
                ## OK   geneTypes47 <- list("mge","nop","other")
                ## OK   geneTypes48 <- list("mge","other")
            ## CASE F:  mge only
                ## OK   geneTypes11 <- list("abc","jkl","mge","nop","qrs")
                ## OK   geneTypes26 <- list("abc","jkl","mge")
                ## OK   geneTypes44 <-  list("mge","nop","qrs")
                ## OK geneTypes45 <-  list("mge")
            ## CASE G: other only
                ## OK   geneTypes12 <- list("abc","jkl","nop","other","qrs")
                ## OK   geneTypes27 <- list("abc","jkl","nop","other")
                ## OK   geneTypes49 <- list("other","qrs")
                ## OK geneTypes50 <- list("other")
            ## CASE H: all annotated
                ## geneTypes28 <- list("abc","jkl","nop","qrs")
        if (length(notMe)>=1) {
            ## if we have hypothetical proteins
            if (notMe == 1) {
                tempColors <- "#FFFFFF"
            } else  {
                tempColors <- geneColors[1:(notMe-1)]
                tempColors <- append(tempColors, "#FFFFFF")
            }
            ## ok if we have mge & other as well
            if ((length(notMeEither)>=1) && (length(norMe)>=1)) {
                print("case A")
                ## if the mge is right after hyp vs. if it's further down
                if (notMe == (notMeEither-1))  {
                    tempColors <- append(tempColors, "#888888")
                } else {
                    temp2Colors <- geneColors[(notMe):(notMeEither-2)]
                    tempColors <- append(tempColors, temp2Colors)
                    tempColors <- append(tempColors, "#888888")
                }
                ## if "other" is right after MGE, vs. if its further down
                if (notMeEither == (norMe-1)) {
                    tempColors <- append(tempColors, "#DEDEDE")
                } else {
                    temp3Colors <- geneColors[(notMeEither-1):(norMe-3)]
                    tempColors <- append(tempColors, temp3Colors)            
                    tempColors <- append(tempColors, "#DEDEDE")
                }
                if (geneColors[(norMe-3)] != geneColors[length(geneColors)] && length(geneColors > 1)) {
                    temp4Colors <-geneColors[(norMe-2):length(geneColors)]
                    tempColors <- append(tempColors, temp4Colors)
                } else if (length(geneColors) == 1 && notMe == 1 && length(geneColors)>0) {
                    ## if we hit this point, we already assigned hyp, mge, and other, and we didn't assign a color, let's add it
                    if (length(tempColors) < 4) {
                        tempColors <- append(tempColors, geneColors)  
                    }  
                }
            } else if (length(notMeEither)>=1) {
                        print("case B")
                ## if there is an MGE (and a hyp) but not an other
                ## if it's right after hyp vs. if it's further down
                if (notMe == (notMeEither-1)) {
                    tempColors <- append(tempColors, "#888888")
                } else {
                    temp5Colors <- geneColors[(notMe):(notMeEither-2)]
                    tempColors <- append(tempColors, temp5Colors)            
                    tempColors <- append(tempColors, "#888888")
                }
                if (geneColors[(notMeEither-2)] != geneColors[length(geneColors)] && length(geneColors) > 1) {
                    temp6Colors <-geneColors[(notMeEither-1):length(geneColors)]
                    tempColors <- append(tempColors, temp6Colors)
                } else if (length(geneColors) == 1 && notMe == 1 && length(geneColors)>0) {
                    ## if we hit this point, we assigned hyp and mge and not the one color
                    if (length(tempColors) < 3) {
                        tempColors <- append(tempColors, geneColors)  
                    }    
                }             
            } else if (length(norMe)>=1) {
                        print("case C")
                ## if there is an MGE (and a hyp) but not an other
                ## if it's right after hyp vs. if it's further down
                if (notMe == (norMe-1)) {
                    tempColors <- append(tempColors, "#DEDEDE")
                } else {
                    temp7Colors <- geneColors[(notMe):(norMe-2)]
                    tempColors <- append(tempColors, temp7Colors)            
                    tempColors <- append(tempColors, "#DEDEDE")
                }
                if (geneColors[(norMe-2)] != geneColors[length(geneColors)] && length(geneColors) > 1) {
                    temp8Colors <-geneColors[(norMe-1):length(geneColors)]
                    tempColors <- append(tempColors, temp8Colors)
                } else if (length(geneColors) == 1 && notMe == 1 && length(geneColors)>0) {
                    ## if we hit this point, we assigned hyp and other and not the one color
                    if (length(tempColors) < 3) {
                        tempColors <- append(tempColors, geneColors)  
                    }  
                }                 
            } else {
                print("case D")
                        ## we've only got hypotheticals
                if (notMe == 1) {
                    ## we can only have assigned hyp so far
                    if (length(geneColors)>0) {
                        tempColors <- append(tempColors, geneColors)
                    }
                } else {
                ## we have already assigned a color + hyp
                    if (length(geneColors)>1) {
                        temp9Colors <-geneColors[(notMe):length(geneColors)]
                        tempColors <- append(tempColors, temp9Colors)
                    }
                }    
            }
        } else {
            ## if we have entered a brave new world where there are no hypothetical proteins
            ## if we have both mge + other
            if ((length(notMeEither)>=1 && length(norMe)>=1)) {
                print("case E")
                ## if MGE is first
                if (notMeEither == 1) {
                    tempColors <- "#888888"
                } else  {
                    tempColors <- geneColors[1:(notMeEither-1)]
                    tempColors <- append(tempColors, "#888888")
                }
                ## if MGE is immediately followed by other
                if (notMeEither == (norMe-1)) {
                    tempColors <- append(tempColors, "#DEDEDE")
                } else {
                    temp2Colors <- geneColors[(notMeEither):(norMe-2)]
                    tempColors <- append(tempColors, temp2Colors)            
                    tempColors <- append(tempColors, "#DEDEDE")
                }
                if (geneColors[(norMe-2)] != geneColors[length(geneColors)] && length(geneColors)>1) {
                    temp3Colors <- geneColors[(norMe-1):length(geneColors)]
                    tempColors <- append(tempColors, temp3Colors)
                } else if (length(geneColors) == 1 && notMeEither == 1) {
                    ## if we hit this point, we assigned mge and other and not the one color
                    if (length(tempColors) < 3) {
                        tempColors <- append(tempColors, geneColors)  
                    }    
                }   
            } else if (length(notMeEither)>=1) {
                print("case F")
                ## if there is an MGE but no other
                if (notMeEither == 1) {
                    tempColors <- "#888888"
                    tempColors <- append(tempColors, geneColors)
                } else  {
                    tempColors <- geneColors[1:(notMeEither-1)]
                    tempColors <- append(tempColors, "#888888")
                    if (geneColors[(notMeEither-1)] != geneColors[length(geneColors)] && length(geneColors)>1) {
                        temp3Colors <- geneColors[(notMeEither):length(geneColors)]
                        tempColors <- append(tempColors, temp3Colors)
                    } else if (length(geneColors) == 1 && notMeEither == 1) {
                    ## if we hit this point, we assigned mge and other and not the one color
                        if (length(tempColors) < 2) {
                            tempColors <- append(tempColors, geneColors)  
                        }    
                    }   
                }              
            } else if (length(norMe)>=1) {
                print("case G")
                ## if there is not an MGE or hyp, but there is still "other"
                if (norMe == 1) {
                    tempColors <- "#DEDEDE"
                    tempColors <- append(tempColors, geneColors)
                } else  {
                    tempColors <- geneColors[1:(norMe-1)]
                    tempColors <- append(tempColors, "#DEDEDE")
                    if (geneColors[(norMe-1)] != geneColors[length(geneColors)] && length(geneColors)>1) {
                        temp3Colors <- geneColors[(norMe):length(geneColors)]
                        tempColors <- append(tempColors, temp3Colors)
                    } else if (length(geneColors) == 1 && norMe == 1) {
                    ## if we hit this point, we assigned mge and other and not the one color
                        if (length(tempColors) < 2) {
                            tempColors <- append(tempColors, geneColors)  
                        }    
                    }   
                }              
            } else {
                print("case H")
                ## if literally everything is annotated somehow
                tempColors <- geneColors
            }
        }
        finalColors <- tempColors
        print("Palettes generated.")
		return(list(finalColors, allGeneTypes))
	}
