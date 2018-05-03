
get_chromosome_length<-function(chr) {
  if(chr == "1") { return (249250621) }
  if(chr == "2") { return (243199373) }
  if(chr == "3") { return (198022430) }
  if(chr == "4") { return (191154276) }
  if(chr == "5") { return (180915260) }
  if(chr == "6") { return (171115067) }
  if(chr == "7") { return (159138663) }
  if(chr == "8") { return (146364022) }
  if(chr == "9") { return (141213431) }
  if(chr == "10") { return (135534747) }
  if(chr == "11") { return (135006516) }
  if(chr == "12") { return (133851895) }
  if(chr == "13") { return (115169878) }
  if(chr == "14") { return (107349540) }
  if(chr == "15") { return (102531392) }
  if(chr == "16") { return (90354753) }
  if(chr == "17") { return (83257441) }
  if(chr == "18") { return (78077248) }
  if(chr == "19") { return (59128983) }
  if(chr == "20") { return (63025520) }
  if(chr == "21") { return (48129895) }
  if(chr == "22") { return (51304566) }
  if(chr == "X") { return (155270560) }
  if(chr == "Y") { return (59373566) }

  return (0)
}

get_centromere_position<-function(chr, arm) {

  if(chr == "1") { if(arm == 'P') { return (121535434) } else { return (124535434) } }
  if(chr == "2") { if(arm == 'P') { return (92326171) } else { return (95326171) } }
  if(chr == "3") { if(arm == 'P') { return (90504854) } else { return (93504854) } }
  if(chr == "4") { if(arm == 'P') { return (49660117) } else { return (52660117) } }
  if(chr == "5") { if(arm == 'P') { return (46405641) } else { return (49405641) } }
  if(chr == "6") { if(arm == 'P') { return (58830166) } else { return (61830166) } }
  if(chr == "7") { if(arm == 'P') { return (58054331) } else { return (61054331) } }
  if(chr == "8") { if(arm == 'P') { return (43838887) } else { return (46838887) } }
  if(chr == "9") { if(arm == 'P') { return (47367679) } else { return (50367679) } }
  if(chr == "10") { if(arm == 'P') { return (39254935) } else { return (42254935) } }
  if(chr == "11") { if(arm == 'P') { return (51644205) } else { return (54644205) } }
  if(chr == "12") { if(arm == 'P') { return (34856694) } else { return (37856694) } }
  if(chr == "13") { if(arm == 'P') { return (16000000) } else { return (19000000) } }
  if(chr == "14") { if(arm == 'P') { return (16000000) } else { return (19000000) } }
  if(chr == "15") { if(arm == 'P') { return (17000000) } else { return (20000000) } }
  if(chr == "16") { if(arm == 'P') { return (35335801) } else { return (38335801) } }
  if(chr == "17") { if(arm == 'P') { return (22263006) } else { return (25263006) } }
  if(chr == "18") { if(arm == 'P') { return (15460898) } else { return (18460898) } }
  if(chr == "19") { if(arm == 'P') { return (24681782) } else { return (27681782) } }
  if(chr == "20") { if(arm == 'P') { return (26369569) } else { return (29369569) } }
  if(chr == "21") { if(arm == 'P') { return (11288129) } else { return (14288129) } }
  if(chr == "22") { if(arm == 'P') { return (13000000) } else { return (16000000) } }
  if(chr == "X") { if(arm == 'P') { return (58632012) } else { return (61632012) } }
  if(chr == "Y") { if(arm == 'P') { return (10104553) } else { return (13104553) } }

  return (0)
}

get_arm_length<-function(chr, arm) {

  centromerePosition = get_centromere_position(chr, arm)

  if(arm == 'P')
  {
    return (centromerePosition)
  }
  else
  {
    chrLength = get_chromosome_length(chr)
    return (chrLength - centromerePosition)
  }
}

get_poisson_prob<-function(observed,  expected, armCount = 44) {

  # calc cummulative probabilty of seeing 1 less, and take reverse prob of this to
  # get prob of observed value or higher
  cumProb = 1- ppois(observed - 1, expected)

  # scale by number of arms
  probAcrossArms = 1 - ((1 - cumProb)^armCount)
  return (probAcrossArms)
}


