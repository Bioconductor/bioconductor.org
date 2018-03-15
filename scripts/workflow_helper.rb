require 'net/http'
require 'uri'

def get_dcfs(repo, version)	
  url = URI.parse("http://master.bioconductor.org/packages/#{version}/#{repo}/VIEWS")
  req = Net::HTTP::Get.new(url.path)
  # FIXME make sure that the request was successful (returns
  # http status 200)
  res = Net::HTTP.start(url.host, url.port) {|http|
    http.request(req)
  }
  views = res.body
  views.force_encoding('UTF-8')
  view_dcfs = []
  view_lines = views.split("\n")
  dcf = ""
  view_lines.each_with_index do |line, i|
    line = line.force_encoding("UTF-8") if line.respond_to? :force_encoding
    if i == (view_lines.length() -1)
      if !line.empty?
        dcf = dcf + "\n" + line
        line = ""
      end
    end
    if line.empty?
      pdcf = DcfHelper.parse(dcf)
      view_dcfs.push pdcf
      dcf = ""
    else
      if dcf == ""
        dcf = line
      else
        dcf = dcf + "\n" + line
      end
    end
  end
  return clean_dcfs(view_dcfs)
end


def clean_dcfs(dcfs)
  ret = {}
  plural_fields = ["Depends", "Suggests", "Imports", "Enhances", "biocViews",
    "LinkingTo",  "vignettes", "vignetteTitles", "Rfiles", "htmlDocs",
    "htmlTitles"]
  for dcf in dcfs
    dcf = dcf.first if dcf.is_a? Array
    for key in dcf.keys
      if plural_fields.include? key
        dcf[key] = DcfHelper.get_value_as_array(dcf[key])
      end
    end
    key = dcf["Package"]
    ret[key] = dcf
  end
  ret
end

