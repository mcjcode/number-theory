import sys
import glob

BG_RED = '\x1B[41m'
RESET = '\x1b[0m'

def checkfile(fn, print_output=False):
    
    out = ''
    fail = False
    
    with open(fn) as fp:
        for line in fp.readlines():
            if line.rstrip() != line[:-1]:
                fail = True
                out += line.rstrip() + BG_RED + ' '*(len(line[:-1])-len(line.rstrip())) + RESET + '\n'
            else:
                out += line
    
    if fail and print_output:
        print(out)
    return not fail
    
def test_files():
    for fn in glob.glob('*.py'):
        assert checkfile(fn)
        
if __name__=='__main__':
    
    fn = sys.argv[1]
    checkfile(fn, print_output=True)