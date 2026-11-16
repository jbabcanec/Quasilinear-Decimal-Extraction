"""Step 1b: same test on the independent v101 transliteration."""
from v101parse import v101_table
from run1_base import run

d = v101_table()
print('v101 numeric-locus lines', len(d), 'paragraphs', d.pid.nunique())
body = d[d.li >= 2]
print('non-first lines', len(body))
print('\n== v101: first line excluded, first v101 char, within-paragraph shuffle ==')
run(body, 'g', label='v101, excl line1')
print('\n== v101 first line included ==')
run(d, 'g', lags=(1, 2), label='v101, incl line1')
print('\n== v101 negative controls ==')
run(body, 'w2g', lags=(1, 2), label='v101 2nd word first glyph')
run(body, 'lastg', lags=(1, 2), label='v101 line-last glyph')
