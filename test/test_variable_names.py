import re
import glob


def test_no_doubled_chars():
    pattern = re.compile('\\b([_a-zA-Z])\\1\\b')

    for fn in glob.glob('../number-theory/*.py'):
        with open(fn) as fp:
            text = fp.read()
            match = pattern.search(text)
            if match:
                print(fn, match)
            #assert match==None, '%s %s'%(fn, match)

if __name__=='__main__':
    test_no_doubled_chars()
