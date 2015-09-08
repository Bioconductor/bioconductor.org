require 'kramdown'

class ImageFilter < Nanoc::Filter
  identifier :image_captions
  
  def run(content, params={})
    content.gsub(/!\[([^\]]+)\]\(([^\)]+)\)/m) do
      alt = $1
      src = $2
      fig = "![#{alt}](#{src})"
     
      fig << "\n<div class=\"caption\">#{Kramdown::Document.new(alt).to_html}</div>" unless alt.empty?
      fig

    end
  end
  
end
