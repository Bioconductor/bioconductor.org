# encoding: utf-8

require 'rubygems'
require 'httparty'
require 'time'
require 'yaml'
require 'fileutils'
require 'nokogiri'
require 'open-uri'
require 'date'

class PipermailList < Nanoc3::DataSource
    identifier :pipermail_list

    def fetch
        baseurl = "https://stat.ethz.ch/pipermail/bioconductor"
        this_month = Date.today.strftime("%B")
        this_year = Date.today.strftime("%Y")
        last_month = (Date.today<<1).strftime("%B")
        year_last_month = (Date.today<<1).strftime("%Y")
        url = "#{baseurl}/#{this_year}-#{this_month}/date.html"
        # TODO ensure url exists?
        # why is the certificate not valid from dan's home laptop?
        doc = Nokogiri::HTML(open(url, {ssl_verify_mode: OpenSSL::SSL::VERIFY_NONE}))
        uls = doc.css("ul")
        return [] if uls.empty?
        msgul = uls[1]
        return [] if msgul.empty?
        a = msgul.css("a")
        return [] if a.empty?
        num_wanted = 5
        index_of_first = nil
        if a.length > num_wanted * 2
            index_of_first = a.length - (num_wanted * 2) - 1
        else
            index_of_first = 0
        end
        ret = []
        for i in index_of_first..a.length-1
            next unless i % 2 == 0
            x = a[i]
            href = x.attr('href')
            subject = x.text.strip.sub(/^\[BioC\] /, "")
            msg_url = url.sub("date.html", href)
            msg_doc = Nokogiri::HTML(open(msg_url, {ssl_verify_mode: OpenSSL::SSL::VERIFY_NONE}))
            datestamp = msg_doc.css("i").first.text#.gsub(/[ ]{2,}/," ")
            dt = DateTime.strptime(datestamp, "%a %b %e %H:%M:%S %Z %Y")
            isodate = dt.iso8601
            puts i
            puts a[i]
        end
    end


    def items
        fetch()
    end

end
