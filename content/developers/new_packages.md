This is a list of the last 100 packages added to Bioconductor and available in the
[development version](/packages/devel/bioc) of Bioconductor. The list is also available
as an [RSS Feed](/rss/new_packages.rss).

<ul>
  <% for package in recent_packages() %>
    <li><a class="symlink" href="/packages/devel/bioc/html/<%=package%>.html"><%=package%></a></li>
  <% end %>
</ul>
