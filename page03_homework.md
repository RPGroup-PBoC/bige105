---
layout: page
title: Homework
img: code.png # Add image post (optional)
permalink: homework
sidebar: true
---

---


<table>
<tr>
<th> <b>Homework</b></th>
<th> <b> Problem Set </b></th>
<th> <b>Associated Files</b></th>
<th> <b> Due Date</b> </th>
</tr>
{% for hwk in site.data.homework %}
<tr>
<td> {{hwk.number}} </td>
<td> <a href="http://www.rpgroup.caltech.edu/bige105/hwk/{{hwk.pset}}"> Problem Set </a></td>
{% if hwk.data %}
{% for dataset in hwk.data %}
<td> <a href="http://www.rpgroup.caltech.edu/bige105/data/{{dataset.link}}">{{dataset.name}}</a></td>
{% endfor %}
{% else %}
<td> -- </td>
{% endif %}
<td> {{hwk.due_date}} </td>
{% if hwk.solutions %}
<td> <a href="https://rpdata.caltech.edu/courses/bige105/2020/{{hwk.solns}}">Solutions</a></td>
{% endif %}
</tr>
{%endfor%}
</table>
