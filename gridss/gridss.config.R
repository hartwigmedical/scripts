# configuration
gridss.short_event_size_threshold = 1000

# somatic filters
gridss.allowable_normal_contamination=0.005 # 0.5% of the supporting fragments can come from the normal
gridss.min_normal_depth = 8

# initial consideration filters
gridss.min_breakpoint_depth = 1
gridss.max_homology_length = 50
gridss.max_allowable_shot_event_strand_bias = 0.95

# final output filters
gridss.min_breakpoint_qual = 500
gridss.min_breakend_qual = 1000
