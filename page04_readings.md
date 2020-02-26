---
layout: page
title: Readings
img: readings.png # Add image post (optional)
permalink: readings
sidebar: true
---

---

As the course progresses, we will post papers we find interesting or that we
specifically mention in class.

{% for pub in site.data.pubs %}
<article class="post">
<a class="post-thumbnail" style="background-image: url(http://rpgroup.caltech.edu/bige105/assets/img/{{pub.pic}})" href="http://rpgroup.caltech.edu/bige105/assets/papers/{{pub.file}}"> </a>
<div class="post-content">
<b class="post-title"><a href="http://rpgroup.caltech.edu/bige105/assets/papers/{{pub.file}}">{{pub.title}}</a></b>
<p>by {{pub.authors}} in <i>{{pub.journal}}</i>, {{pub.vol_iss}} {{pub.year}}.</p>
<p>•<a href="{{pub.publisher_link}}">Publisher</a><br/></p>
<p>•<a href="http://rpgroup.caltech.edu/bige105/assets/papers/{{pub.file}}">PDF</a><br/></p>
</div>
</article>
{%endfor%}
