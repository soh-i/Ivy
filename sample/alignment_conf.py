from Ivy.alignment import AlignmentConfig

conf = AlignmentConfig()
conf.print_all_params()

#print conf.conf
print conf.get_filter_value("mapq")
print conf.has_filter('is_duplicate')
print conf.set_filter('base_qual', 33)
print conf.get_filter_value("base_qual")
