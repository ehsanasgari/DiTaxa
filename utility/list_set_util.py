
import operator

def get_intersection_of_list(list_of_list_features):
    return list(set.intersection(*map(set, list_of_list_features)))

def get_max_of_dict(inp):
    return max(inp.items(), key=operator.itemgetter(1))[0]
