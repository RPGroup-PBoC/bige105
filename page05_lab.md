---
layout: page
title: Laboratory Experiments
img: lab.png # Add image post (optional)
permalink: lab
sidebar: true
---

This course will include a laboratory section where we will observe evolution in real time.

{% for lab in site.data.labs %}
<article class="post">
<a class="post-thumbnail" style="background-image: url(http://rpgroup.caltech.edu/bige105/assets/img/{{lab.pic}})" href="{{site.url}}/{{site.baseurl}}/assets/labs/{{lab.file}}"> </a>
<div class="post-content">
<b class="post-title"><a href="http://rpgroup.caltech.edu/bige105/assets/labs/{{lab.file}}">{{lab.title}}</a></b>
<p>•<a href="http://rpgroup.caltech.edu/bige105/assets/labs/{{lab.file}}">PDF</a><br/></p>
</div>
</article>
{%endfor%}
