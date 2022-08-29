import re

invalid = re.finditer("[^ATCG]+", "TGFA")
matches = len([invalid])

print(i
print(matches)