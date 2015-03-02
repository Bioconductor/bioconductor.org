![](/images/icons/magnifier.gif)Mirrors
=======================================

<% for country in config[:mirrors] %>
<%= country.keys.first.to_s %>
------------------------

<% for mirror in country.values.first %>
* [<%= mirror[:institution] %>](<%= mirror[:institution_url] %>)

  URL: <<%= mirror[:mirror_url] %>>

  Contact: <%= render_mirror_contacts(mirror) %>

<% end %>
<% end %>

