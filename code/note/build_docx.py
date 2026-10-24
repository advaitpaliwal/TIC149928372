# Build rnaas11.docx: each author-block line and each reference becomes its own paragraph (pandoc joins adjacent lines).
import re, pypandoc
s = open('rnaas11.md').read()
head, rest = s.split('## Abstract', 1)
body, refs = rest.split('## References', 1)
head = re.sub(r'\n(?=\S)', '\n\n', head.strip()) + '\n\n'
refs = '\n\n'.join(l for l in refs.strip().splitlines() if l.strip())
open('glue/_docx_src.md', 'w').write(head + '## Abstract' + body + '## References\n\n' + refs + '\n')
pypandoc.convert_file('glue/_docx_src.md', 'docx', outputfile='rnaas11.docx', extra_args=['--resource-path=.'])
