Mirrors
=======================================

If you are interested in maintaining a mirror of this site (for either
public or private use) [read this](mirror-how-to/).

<% for country in config[:mirrors] %>
  <% next if country.keys.first.to_s == "0-Bioconductor" %>
<%= country.keys.first.to_s %>
------------------------

<% for mirror in country.values.first %>
* [<%= mirror[:institution] %>](<%= mirror[:institution_url] %>)

  URLs: <%= render_mirror_urls(mirror) %>

  Contact: <%= render_mirror_contacts(mirror) %>

<% end %>
<% end %>
