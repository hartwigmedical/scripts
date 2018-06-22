# configuration
gridss.short_event_size_threshold = 1000

# somatic filters
gridss.allowable_normal_contamination=0.005 # 0.5% of the supporting fragments can come from the normal
gridss.min_normal_depth = 8

# initial consideration filters
gridss.min_direct_read_support = 5
gridss.max_homology_length = 50
gridss.max_allowable_short_event_strand_bias = 0.95
gridss.single_breakend_multiplier = 2 # Require this much more support for single breakends

# final output filters
gridss.min_qual = 500
