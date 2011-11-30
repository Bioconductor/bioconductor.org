from HTMLParser import HTMLParser
import time
import sys
import subprocess

global data_dict
data_dict = {}


class MyHTMLParser(HTMLParser):
    want_data = False
    
    def handle_starttag(self, tag, attrs):
        name = None
        checked = False
        for attr in attrs:
            if attr[0] == 'name':
                name = attr[1]
            if attr[0] == 'checked':
                checked = True
        
        if tag == 'textarea' or checked:
            self.want_data = name
            
    def handle_data(self, data):
        if (self.want_data != False):
            data_dict[self.want_data] = data
            self.want_data = False

if __name__ == "__main__":
    argc = len(sys.argv)
    if (argc != 2):
        print "supply an email address."
        sys.exit(1)
    
    email = sys.argv[1]
    
    ts = time.localtime()
    secs = int(time.mktime(ts))

    cookiefile = "/tmp/upl_cookies_%d.txt" % secs
    
    parser = MyHTMLParser()
    
    url = "https://stat.ethz.ch/mailman/admin/bioconductor"
    #url = "http://localhost/~dtenenba/test.cgi"
    
    pwf = open("wordfile")
    pw = pwf.read().strip()
    pwf.close
    
    cmd = 'curl -s -o /dev/null --cookie-jar %s  -F "adminpw=%s" %s'\
      % (cookiefile, pw, url)
    #print cmd
    retcode = subprocess.call(cmd, shell=True)
    if (retcode != 0):
        print "failed to log in"
        sys.exit(retcode)
    
    htmlfile = "/tmp/upl_%d_forminfo.html" % secs
    cmd = "curl -s --cookie %s -o %s https://stat.ethz.ch/mailman/admin/bioconductor/privacy/sender" \
      % (cookiefile, htmlfile)
    #print cmd
    retcode = subprocess.call(cmd, shell=True)
    if (retcode != 0):
        print "failed to get form info"
        sys.exit(retcode)
     
    f = open(htmlfile)
    html = f.read()
    f.close()

    parser.feed(html)

    #print "data dict is:"
    #print data_dict

    
    cmd = "curl -s --cookie %s " % cookiefile

    for key, value in data_dict.items():
        filename = "/tmp/upl_%d_%s.txt" % (secs, key)
        f = open(filename, 'w')
        if key == 'accept_these_nonmembers':
            value += "\n%s" % email
            f.write(value)
        f.close
        cmd += '-F "%s=<%s" ' % (key, filename)

    #cmd += " http://localhost/~dtenenba/test.cgi"
    cmd += " https://stat.ethz.ch/mailman/admin/bioconductor/privacy/sender"

    #print cmd
    
    ## TODO - cleanup /tmp/upl* files when done
    
    retcode = subprocess.call(cmd, shell=True)
    
    ## TODO - logout?
    print "command completed with return code: %d" % retcode
    sys.exit(retcode)
