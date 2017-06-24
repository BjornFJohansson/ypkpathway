from IPython.nbformat import current as nbf

nb = nbf.new_notebook()

cells = []

for var in [1,2,3]:
    text = "Variable{}".format(var)
    cell = nbf.new_text_cell('markdown', text)
    cells.append(cell)
    cell = nbf.new_code_cell("print 123")
    cells.append(cell)

nb['worksheets'].append(nbf.new_worksheet(cells=cells))

with open('generated_notebook.ipynb', 'w') as f:
        nbf.write(nb, f, 'ipynb')
