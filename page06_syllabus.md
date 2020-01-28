---
layout: page
title: Syllabus
img: lab.png # Add image post (optional)
permalink: syllabus
sidebar: true
---

<table>
<tr>
    <th><b>Date</b></th>
    <th><b>Week</b></th>
    <th><b>Topic</b></th>
    <th><b>Due</b></th>
    <th><b>Slides</b></th>
    <th><b>Reading</b></th>
</tr>
{% for day in site.data.syllabus %}
<tr>
    <td>{{day.date}}</td>
    <td>{{day.week}}</td>
    <td>{{day.topic}}</td>
    {% if day.due %}
    <td>{{day.due}}</td>
    {% else %}
    <td> -- </td>
    {% endif %}
    {% if day.slides %}
    <td><a href="http://rpdata.caltech.edu/courses/bige105/{{day.slides}}">
    PDF </a></td>
    {% else %}
    <td> -- </td>
    {% endif %}
    <td>{{day.reading}}</td>
</tr>
{% endfor %}
</table>
