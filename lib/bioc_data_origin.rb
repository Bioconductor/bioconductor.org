require 'yaml'

# Origin for build/check data that the site build CONSUMES.
#
# This is deliberately separate from the URLs the build WRITES into pages. Those
# must stay bioconductor.org — they are links for readers. This is the other
# direction: where the build reads VIEWS, BUILD_STATUS_DB.txt, checkResults
# summaries and the git logs from.
#
# Historically every one of those reads was hardcoded to bioconductor.org or
# master.bioconductor.org, meaning the build fetched its inputs from the site it
# produces. That works only while that server exists and already holds current
# data; it is not a real origin, and it makes the build impossible to run
# anywhere else. Routing the reads through one setting makes the dependency
# visible and repointable without touching code.
#
# Resolution order: BIOC_BUILD_DATA_ORIGIN env var, then build_data_origin in
# config.yaml, then the historical default so behaviour is unchanged by default.
module BiocDataOrigin
  DEFAULT = 'https://master.bioconductor.org'.freeze

  def self.url
    @url ||=
      ENV['BIOC_BUILD_DATA_ORIGIN'] ||
      (File.exist?('./config.yaml') && YAML.load_file('./config.yaml')['build_data_origin']) ||
      DEFAULT
  end
end

def bioc_data_origin
  BiocDataOrigin.url
end
