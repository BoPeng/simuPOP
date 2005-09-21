a="(int a, int b=vector(), double c=3)"
b=a.split(',')
out=[]
for s in b:
  piece = s.split('=')
  var = piece[0].split(' ')[-1]
  if( len( piece ) == 2 ):   # has default
      defVal = piece[1].split('(')[0].split(')')[0]
      defVal = defVal.replace('vector','[]')
      out.append( var + '=' + defVal  )
  else:
      out.append( var )

print '(' + (', '.join(out)) + ')'
  
help(population.__init__)
a=5 
