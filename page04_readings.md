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
<li> <a style="text-decoration: none;" href="https://rpdata.caltech.edu/2020/protected/papers/{{pub.fname}}"> <b><i>{{pub.title}}</i> by {{pub.authors}} in {{pub.journal}}, {{pub.year}}, {{pub.vol_iss}}.</b></a> {{pub.desc}}</li>
{% endfor %}