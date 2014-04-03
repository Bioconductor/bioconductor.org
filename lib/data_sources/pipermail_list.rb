# encoding: utf-8

require 'rubygems'
require 'httparty'
require 'time'
require 'yaml'
require 'fileutils'
require 'nokogiri'
require 'openuri'

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
        doc = Nokogiri::HTML(open(url))
        uls = doc.css("ul")
        return [] if uls.empty?
        msgul = uls[1]
        return [] if msgul.empty?
        a = msgul.css("a")
        return [] if a.empty?
        num_wanted = 5
        

    end


    def items
        fetch()
    end

end
