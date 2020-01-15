---
layout: page
title: Computational Explorations
img: code.png # Add image post (optional)
permalink: code
sidebar: true
---

---

During this course, you will develop a computational prowess that will aid in
your understanding of evolution. We will post Jupyter Notebooks of the tutorial
sessions here. 


{% if site.data.tutorials %}
{% for tut in site.data.tutorials %}
<article class="post">
<a class="post-thumbnail" style="background-image: url(http://rpgroup.caltech.edu/bige105/assets/img/{{tut.pic}})" href="{{site.url}}/{{site.baseurl}}/tutorials/{{tut.link}}.html"> </a>

<div class="post-content">
<b class="post-title"><a href="http://rpgroup.caltech.edu/bige105/tutorials/{{tut.link}}.html">{{tut.title}}</a></b>
<p> {{tut.desc}}</p>
<p>• <a href="http://rpgroup.caltech.edu/bige105/tutorials/{{tut.link}}.ipynb"> Jupyter Notebook (.ipynb)</a><br/></p>
{% if tut.req %}<i>Necessary Data Sets </i><br/>
{% for ds in tut.req %}
{% if ds.storage == 'local' %}
{% assign link = "{{ds.link}}" %}
{% else %}
{% assign link = "{{ds.link}}" %}
{% endif %}
<a style="font-size: 0.9em;" href="{{link}}"> - {{ds.title}} </a><br/>
{% endfor %}
{% endif %}
</div>
</article>
{%endfor%}
{% endif %}