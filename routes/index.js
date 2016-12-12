var express = require('express');
var router = express.Router();

/* GET home page. */
router.get('/', function(req, res, next) {
  res.render('index', { title: 'sb-blog' });
});
router.get('/article',(req,res)=>{
  res.render('/article')
})

module.exports = router;
