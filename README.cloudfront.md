# CloudFront

[CloudFront](https://aws.amazon.com/cloudfront/) is a content delivery
network (CDN). The idea is that if a user in say, Singapore hits our web
site, the files she requested will be cached in a location near her and her next access (or that of any other user in that area) will be faster.

Our CloudFront distribution can be managed through the
[AWS Console](https://console.aws.amazon.com/cloudfront/home?region=us-east-1#).

Once CloudFront is running, there isn't much you have to do.

Just be aware that when you are outside the Hutch network and you 
hit bioconductor.org, you are hitting CloudFront.

A good way to test this is as follows:

```
$ curl -I http://bioconductor.org/
HTTP/1.1 200 OK
Content-Type: text/html
Content-Length: 16531
Connection: keep-alive
Date: Tue, 19 May 2015 17:38:25 GMT
Server: Apache/2.2.12 (Linux/SUSE)
Last-Modified: Tue, 19 May 2015 17:30:47 GMT
ETag: "62f141-4093-51672abbe1fc0"
Accept-Ranges: bytes
Cache-Control: max-age=600
Expires: Tue, 19 May 2015 17:48:25 GMT
Vary: Accept-Encoding
X-Cache: Miss from cloudfront
Via: 1.1 70d79aa19e315b47281005f9e3c25c88.cloudfront.net (CloudFront)
X-Amz-Cf-Id: 97PFtz_9dGUDk2kHEjggRwbhRkcSG0vNg9akW2y70mvyK0iKO882MA==
```

The `X-cache` and `X-Amz...` headers tell you you are dealing with CloudFront.
If you want to hit the web server at FHCRC, go instead to
[http://master.bioconductor.org](http://master.bioconductor.org).

Different file types have different cache expiration times, this
is set in the 
[.htaccess file](https://hedgehog.fhcrc.org/bioconductor/trunk/bioconductor.org/assets/.htaccess) 
for the web site.

If a file gets propagated by CloudFront that shouldn't be, you can remove
it by invalidating it, go to 
[the console page for our distribution](https://console.aws.amazon.com/cloudfront/home?region=us-east-1#distribution-settings:E1TVLJONPTUXV3) and click on the "Invalidations" tab.

