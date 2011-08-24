This is a list of the last 100 packages added to Bioconductor and available in the
[development version](/packages/devel/bioc) of Bioconductor. The list is also available
as an [RSS Feed](/rss/new_packages.rss).

<% for package in recent_packages() %>
  <a class="symlink" href="/packages/devel/bioc/html/<%=package[:package]%>.html"><%=package[:package]%></a> <%= package[:title]%>

><%= package[:description]%>
<% end %>
