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
<th> <b>Description</b></th>
<th> <b> Due Date</b> </th>
<th> <b> Solutions</b> </th><br/>
</tr>
{% for hwk in site.data.homework %}
<tr>
<td> <a href="hwk/{{hwk.pset}}">Problem Set </a></td>
<td> {{hwk.desc}} </td> 
<td> {{hwk.due_date}} </td> 
<td> <a href="https://rpdata.caltech.edu/courses/bige105/2020/{{hwk.solns}}">Solutions</a></td>
</tr>
{%endfor%}
</table>


